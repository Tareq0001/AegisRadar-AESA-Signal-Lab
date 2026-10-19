/**
 * ============================================================================
 * AegisRadar AESA DSP Engine — High-Performance Fast Fourier Transform (FFT)
 * Linear Frequency Modulation (LFM) Chirp Synthesis & Matched Filter Pulse Compression
 * ============================================================================
 * Architected for real-time coherent Pulse-Doppler radar signal processing.
 */

export class Complex {
  constructor(public re: number = 0, public im: number = 0) {}

  static fromPolar(r: number, theta: number): Complex {
    return new Complex(r * Math.cos(theta), r * Math.sin(theta));
  }

  static zero(): Complex {
    return new Complex(0, 0);
  }

  static fromReal(re: number): Complex {
    return new Complex(re, 0);
  }

  add(other: Complex): Complex {
    return new Complex(this.re + other.re, this.im + other.im);
  }

  sub(other: Complex): Complex {
    return new Complex(this.re - other.re, this.im - other.im);
  }

  mul(other: Complex): Complex {
    return new Complex(
      this.re * other.re - this.im * other.im,
      this.re * other.im + this.im * other.re
    );
  }

  scale(scalar: number): Complex {
    return new Complex(this.re * scalar, this.im * scalar);
  }

  div(other: Complex): Complex {
    const denom = other.re * other.re + other.im * other.im;
    if (denom === 0) throw new Error('Complex division by zero');
    return new Complex(
      (this.re * other.re + this.im * other.im) / denom,
      (this.im * other.re - this.re * other.im) / denom
    );
  }

  conj(): Complex {
    return new Complex(this.re, -this.im);
  }

  abs(): number {
    return Math.hypot(this.re, this.im);
  }

  absSq(): number {
    return this.re * this.re + this.im * this.im;
  }

  phase(): number {
    return Math.atan2(this.im, this.re);
  }
}

/**
 * Returns the smallest power of 2 greater than or equal to n.
 */
export function nextPowerOfTwo(n: number): number {
  if (n <= 1) return 1;
  let p = 1;
  while (p < n) {
    p <<= 1;
  }
  return p;
}

/**
 * Radix-2 Cooley-Tukey In-Place Decimation-in-Time Fast Fourier Transform
 * @param input Array of Complex numbers (must have length equal to a power of 2)
 * @param inverse If true, computes the Inverse FFT (IFFT) with 1/N scaling
 */
export function fftRadix2(input: Complex[], inverse: boolean = false): Complex[] {
  const n = input.length;
  if ((n & (n - 1)) !== 0) {
    throw new Error(`FFT length must be a power of 2, received ${n}`);
  }

  // Create working array
  const output: Complex[] = new Array(n);
  for (let i = 0; i < n; i++) {
    output[i] = new Complex(input[i].re, input[i].im);
  }

  // Bit-reversal permutation
  let j = 0;
  for (let i = 0; i < n - 1; i++) {
    if (i < j) {
      const temp = output[i];
      output[i] = output[j];
      output[j] = temp;
    }
    let k = n >> 1;
    while (k <= j) {
      j -= k;
      k >>= 1;
    }
    j += k;
  }

  // Butterfly computations
  const sign = inverse ? 1 : -1;
  for (let len = 2; len <= n; len <<= 1) {
    const halfLen = len >> 1;
    const angle = (sign * 2 * Math.PI) / len;
    const wBase = Complex.fromPolar(1, angle);

    for (let i = 0; i < n; i += len) {
      let w = new Complex(1, 0);
      for (let k = 0; k < halfLen; k++) {
        const u = output[i + k];
        const v = output[i + k + halfLen].mul(w);
        output[i + k] = u.add(v);
        output[i + k + halfLen] = u.sub(v);
        w = w.mul(wBase);
      }
    }
  }

  // Scaling for inverse FFT
  if (inverse) {
    for (let i = 0; i < n; i++) {
      output[i] = output[i].scale(1 / n);
    }
  }

  return output;
}

export interface LFMChirpParams {
  pulseWidth: number;      // Tau (seconds), e.g. 10e-6 s
  bandwidth: number;       // B (Hz), e.g. 20e6 Hz
  samplingRate: number;    // Fs (Hz), e.g. 50e6 Hz
  centerFreq?: number;     // Carrier/IF freq (Hz), default 0 (baseband)
}

export interface LFMChirpResult {
  samples: Complex[];
  time: Float64Array;
  numSamples: number;
  timeBandwidthProduct: number;
  rangeResolutionMeters: number;
  uncompressedResolutionMeters: number;
}

const SPEED_OF_LIGHT = 299792458; // m/s

/**
 * Synthesizes an analytic complex baseband Linear Frequency Modulated (LFM) chirp.
 * Waveform: s(t) = exp(j * pi * K * t^2) where K = B / tau.
 */
export function generateLFMChirp(params: LFMChirpParams): LFMChirpResult {
  const { pulseWidth, bandwidth, samplingRate, centerFreq = 0 } = params;
  const numSamples = Math.floor(pulseWidth * samplingRate);
  const dt = 1 / samplingRate;
  const chirpRate = bandwidth / pulseWidth; // Hz/s

  const time = new Float64Array(numSamples);
  const samples: Complex[] = new Array(numSamples);

  for (let i = 0; i < numSamples; i++) {
    // Centered time window t in [-tau/2, +tau/2]
    const t = -pulseWidth / 2 + i * dt;
    time[i] = t;

    // Instantaneous phase phi(t) = 2*pi*fc*t + pi*K*t^2
    const phase = 2 * Math.PI * centerFreq * t + Math.PI * chirpRate * t * t;
    samples[i] = Complex.fromPolar(1.0, phase);
  }

  const timeBandwidthProduct = bandwidth * pulseWidth;
  const rangeResolutionMeters = SPEED_OF_LIGHT / (2 * bandwidth);
  const uncompressedResolutionMeters = (SPEED_OF_LIGHT * pulseWidth) / 2;

  return {
    samples,
    time,
    numSamples,
    timeBandwidthProduct,
    rangeResolutionMeters,
    uncompressedResolutionMeters,
  };
}

/**
 * Matched Filter Pulse Compression via Frequency-Domain Fast Convolution
 * Computes y[n] = IFFT(FFT(rxSignal) * conj(FFT(referencePulse)))
 * Yields matched filter output with pulse compression gain = 10 * log10(Time-Bandwidth)
 */
export function matchedFilterPulseCompression(
  rxSignal: Complex[],
  referencePulse: Complex[]
): {
  compressedSignal: Complex[];
  magnitudes: Float64Array;
  peakIndex: number;
  peakMagnitude: number;
} {
  const minLen = rxSignal.length + referencePulse.length - 1;
  const fftLen = nextPowerOfTwo(minLen);

  // Zero-pad received signal to fftLen
  const paddedRx: Complex[] = new Array(fftLen);
  for (let i = 0; i < fftLen; i++) {
    paddedRx[i] = i < rxSignal.length ? rxSignal[i] : Complex.zero();
  }

  // Zero-pad reference pulse to fftLen
  const paddedRef: Complex[] = new Array(fftLen);
  for (let i = 0; i < fftLen; i++) {
    paddedRef[i] = i < referencePulse.length ? referencePulse[i] : Complex.zero();
  }

  // Forward FFTs
  const rxFreq = fftRadix2(paddedRx, false);
  const refFreq = fftRadix2(paddedRef, false);

  // Matched filter multiplication: H_opt(f) = S*(f)
  const productFreq: Complex[] = new Array(fftLen);
  for (let i = 0; i < fftLen; i++) {
    productFreq[i] = rxFreq[i].mul(refFreq[i].conj());
  }

  // IFFT to reconstruct compressed time-domain pulse
  const compressedSignal = fftRadix2(productFreq, true);

  // Extract magnitudes & locate detection peak
  const outLen = rxSignal.length;
  const magnitudes = new Float64Array(outLen);
  let peakIndex = 0;
  let peakMagnitude = -1;

  for (let i = 0; i < outLen; i++) {
    const mag = compressedSignal[i].abs();
    magnitudes[i] = mag;
    if (mag > peakMagnitude) {
      peakMagnitude = mag;
      peakIndex = i;
    }
  }

  return {
    compressedSignal: compressedSignal.slice(0, outLen),
    magnitudes,
    peakIndex,
    peakMagnitude,
  };
}

/**
 * Computes theoretical Pulse Compression Ratio (PCR) in decibels
 */
export function computePulseCompressionGainDb(timeBandwidthProduct: number): number {
  if (timeBandwidthProduct <= 0) return 0;
  return 10 * Math.log10(timeBandwidthProduct);
}
