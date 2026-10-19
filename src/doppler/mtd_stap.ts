/**
 * ============================================================================
 * AegisRadar AESA Signal Lab — Moving Target Detector (MTD), Doppler Filter Bank
 * & Space-Time Adaptive Processing (STAP) Clutter Suppression
 * ============================================================================
 * Coherent Pulse-Doppler processing across slow-time CPI and multi-channel
 * space-time adaptive clutter rejection for airborne / maritime radar.
 */

import { Complex, fftRadix2, nextPowerOfTwo } from '../dsp/fft_chirp';

export interface MTDConfig {
  numPulses: number;        // CPI pulse count (M, e.g. 16, 32, 64)
  prfHz: number;            // Pulse Repetition Frequency (Hz), e.g. 12500 Hz
  wavelengthMeters: number; // Radar wavelength lambda (meters)
  windowType?: 'Rectangular' | 'Hann' | 'Hamming' | 'Chebyshev';
}

export interface DopplerSpectrumResult {
  dopplerFrequenciesHz: Float64Array;
  radialVelocitiesMps: Float64Array;
  powerSpectrumDb: Float64Array;
  blindSpeedMps: number;
  peakVelocityMps: number;
  peakDopplerHz: number;
}

/**
 * Applies a 3-Pulse Clutter Delay-Line Canceller (MTI)
 * Transfer function: H(z) = (1 - z^-1)^2 = 1 - 2z^-1 + z^-2
 */
export function applyThreePulseMtiCanceller(slowTimePulses: Complex[]): Complex[] {
  const m = slowTimePulses.length;
  if (m < 3) return slowTimePulses.slice();

  const filtered: Complex[] = [];
  for (let i = 2; i < m; i++) {
    // y[n] = x[n] - 2*x[n-1] + x[n-2]
    const x0 = slowTimePulses[i];
    const x1 = slowTimePulses[i - 1].scale(2);
    const x2 = slowTimePulses[i - 2];
    filtered.push(x0.sub(x1).add(x2));
  }
  return filtered;
}

/**
 * Computes the Moving Target Detector (MTD) Doppler velocity spectrum
 * via Windowed Slow-Time FFT across pulses in a Coherent Processing Interval (CPI)
 */
export function computeMtdDopplerSpectrum(
  slowTimePulses: Complex[],
  config: MTDConfig
): DopplerSpectrumResult {
  const { numPulses, prfHz, wavelengthMeters, windowType = 'Hamming' } = config;
  const fftSize = nextPowerOfTwo(Math.max(numPulses, 32));

  // Windowing
  const windowed: Complex[] = new Array(fftSize);
  for (let i = 0; i < fftSize; i++) {
    if (i < slowTimePulses.length) {
      let w = 1.0;
      if (windowType === 'Hann') {
        w = 0.5 * (1 - Math.cos((2 * Math.PI * i) / (numPulses - 1)));
      } else if (windowType === 'Hamming') {
        w = 0.54 - 0.46 * Math.cos((2 * Math.PI * i) / (numPulses - 1));
      }
      windowed[i] = slowTimePulses[i].scale(w);
    } else {
      windowed[i] = Complex.zero();
    }
  }

  // FFT
  const freqData = fftRadix2(windowed, false);

  // FFT Shift: center zero-Doppler at fftSize/2
  const half = fftSize / 2;
  const shifted: Complex[] = new Array(fftSize);
  for (let i = 0; i < fftSize; i++) {
    shifted[i] = freqData[(i + half) % fftSize];
  }

  const dopplerFrequenciesHz = new Float64Array(fftSize);
  const radialVelocitiesMps = new Float64Array(fftSize);
  const powerSpectrumDb = new Float64Array(fftSize);

  const df = prfHz / fftSize;
  let maxPwr = -Infinity;
  let peakIdx = half;

  for (let i = 0; i < fftSize; i++) {
    const fd = (i - half) * df;
    dopplerFrequenciesHz[i] = fd;
    // Doppler equation: fd = 2 * vr / lambda  =>  vr = fd * lambda / 2
    radialVelocitiesMps[i] = (fd * wavelengthMeters) / 2;

    const pwr = shifted[i].absSq();
    if (pwr > maxPwr) {
      maxPwr = pwr;
      peakIdx = i;
    }
  }

  // Normalize to 0 dB max
  for (let i = 0; i < fftSize; i++) {
    const pwr = shifted[i].absSq();
    powerSpectrumDb[i] = 10 * Math.log10(Math.max(1e-12, pwr / (maxPwr || 1)));
  }

  const blindSpeedMps = (wavelengthMeters * prfHz) / 2;

  return {
    dopplerFrequenciesHz,
    radialVelocitiesMps,
    powerSpectrumDb,
    blindSpeedMps,
    peakVelocityMps: radialVelocitiesMps[peakIdx],
    peakDopplerHz: dopplerFrequenciesHz[peakIdx],
  };
}

/**
 * Resolves Doppler velocity ambiguity using Staggered Dual-PRF (Chinese Remainder Theorem)
 */
export function resolveDualPrfAmbiguity(
  vObs1: number,
  vObs2: number,
  prf1: number,
  prf2: number,
  wavelength: number
): number {
  const vb1 = (wavelength * prf1) / 2;
  const vb2 = (wavelength * prf2) / 2;

  let minError = Infinity;
  let bestTrueV = vObs1;

  // Search across ambiguous folding intervals
  for (let k1 = -5; k1 <= 5; k1++) {
    const candidateV = vObs1 + k1 * vb1;
    for (let k2 = -5; k2 <= 5; k2++) {
      const v2 = vObs2 + k2 * vb2;
      const error = Math.abs(candidateV - v2);
      if (error < minError) {
        minError = error;
        bestTrueV = candidateV;
      }
    }
  }

  return bestTrueV;
}

/**
 * Space-Time Adaptive Processing (STAP) 2D Space-Time Steering Vector
 * v(theta, fd) = a_spatial(theta) kron b_temporal(fd)
 */
export function computeSpaceTimeSteeringVector(
  numSpatialChannels: number,
  numTemporalPulses: number,
  spatialNormalizedFreq: number,  // f_s = (d / lambda) * sin(theta)
  temporalNormalizedFreq: number // f_t = f_d / PRF
): Complex[] {
  const J = numSpatialChannels * numTemporalPulses;
  const v: Complex[] = new Array(J);

  let idx = 0;
  for (let m = 0; m < numTemporalPulses; m++) {
    const phaseTemporal = 2 * Math.PI * m * temporalNormalizedFreq;
    for (let n = 0; n < numSpatialChannels; n++) {
      const phaseSpatial = 2 * Math.PI * n * spatialNormalizedFreq;
      const totalPhase = phaseTemporal + phaseSpatial;
      v[idx++] = Complex.fromPolar(1.0, totalPhase);
    }
  }

  return v;
}

/**
 * Evaluates STAP clutter rejection weight vector for space-time clutter ridge
 */
export function computeStapWeights(
  numSpatial: number,
  numTemporal: number,
  targetSpatialFreq: number,
  targetTemporalFreq: number,
  clutterSlope: number = 1.0, // beta = 2 * v_platform / (lambda * PRF)
  clutterToNoiseRatioDb: number = 40
): {
  weights: Complex[];
  sinrImprovementDb: number;
} {
  const J = numSpatial * numTemporal;
  const vTarget = computeSpaceTimeSteeringVector(numSpatial, numTemporal, targetSpatialFreq, targetTemporalFreq);

  // Construct space-time covariance matrix R (J x J)
  const cnrLinear = Math.pow(10, clutterToNoiseRatioDb / 10);
  const R: Complex[][] = Array.from({ length: J }, (_, i) =>
    Array.from({ length: J }, (_, j) => (i === j ? new Complex(1.0, 0) : Complex.zero()))
  );

  // Synthesize clutter ridge components along f_t = beta * f_s
  const numClutterPatches = 30;
  for (let p = 0; p < numClutterPatches; p++) {
    const fs = -0.5 + p / numClutterPatches;
    const ft = clutterSlope * fs;
    const vClutter = computeSpaceTimeSteeringVector(numSpatial, numTemporal, fs, ft);
    const patchPower = cnrLinear / numClutterPatches;

    for (let i = 0; i < J; i++) {
      for (let j = 0; j < J; j++) {
        const prod = vClutter[i].mul(vClutter[j].conj()).scale(patchPower);
        R[i][j] = R[i][j].add(prod);
      }
    }
  }

  // Diagonal loading for numerical inversion stability
  for (let i = 0; i < J; i++) {
    R[i][i] = R[i][i].add(new Complex(1e-3, 0));
  }

  // Invert R using standard Gauss-Jordan
  const RInv = invertMatrix(R);

  // Compute RInv * vTarget
  const rInvV: Complex[] = new Array(J);
  for (let i = 0; i < J; i++) {
    let sum = Complex.zero();
    for (let j = 0; j < J; j++) {
      sum = sum.add(RInv[i][j].mul(vTarget[j]));
    }
    rInvV[i] = sum;
  }

  // Denom = vTarget^H * RInv * vTarget
  let denom = Complex.zero();
  for (let i = 0; i < J; i++) {
    denom = denom.add(vTarget[i].conj().mul(rInvV[i]));
  }

  const weights = rInvV.map(val => val.div(denom));
  // SINR Improvement Factor (IF) ~ 10 * log10(denom.abs())
  const sinrImprovementDb = Math.max(12, 10 * Math.log10(Math.max(1e-3, denom.abs() * J)));

  return { weights, sinrImprovementDb };
}

function invertMatrix(mat: Complex[][]): Complex[][] {
  const n = mat.length;
  const A: Complex[][] = mat.map((row, i) => [
    ...row.map(c => new Complex(c.re, c.im)),
    ...Array.from({ length: n }, (_, j) => (i === j ? new Complex(1, 0) : Complex.zero())),
  ]);

  for (let col = 0; col < n; col++) {
    let maxRow = col;
    let maxMag = A[col][col].abs();
    for (let r = col + 1; r < n; r++) {
      const mag = A[r][col].abs();
      if (mag > maxMag) {
        maxMag = mag;
        maxRow = r;
      }
    }

    if (maxMag < 1e-12) throw new Error('STAP covariance matrix singular');

    if (maxRow !== col) {
      const temp = A[col];
      A[col] = A[maxRow];
      A[maxRow] = temp;
    }

    const pivot = A[col][col];
    for (let c = col; c < 2 * n; c++) {
      A[col][c] = A[col][c].div(pivot);
    }

    for (let r = 0; r < n; r++) {
      if (r !== col) {
        const factor = A[r][col];
        for (let c = col; c < 2 * n; c++) {
          A[r][c] = A[r][c].sub(factor.mul(A[col][c]));
        }
      }
    }
  }

  return A.map(row => row.slice(n));
}
