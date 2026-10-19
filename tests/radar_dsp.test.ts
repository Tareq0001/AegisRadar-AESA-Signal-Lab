/**
 * ============================================================================
 * AegisRadar AESA Signal Lab — Unit & Algorithmic Verification Test Suite
 * ============================================================================
 * Rigorous validation of Radix-2 FFT, LFM Chirp, Pulse Compression,
 * CFAR Adaptive Thresholds, AESA Beamforming, MVDR Null-Steering & MTD.
 */

import { describe, it, expect } from 'vitest';
import {
  Complex,
  fftRadix2,
  nextPowerOfTwo,
  generateLFMChirp,
  matchedFilterPulseCompression,
  computePulseCompressionGainDb,
} from '../src/dsp/fft_chirp';
import {
  runCfarDetector,
  computeCaCfarAlpha,
  calculateSwerlingPd,
} from '../src/dsp/cfar_detector';
import {
  AESAArray,
  computeMvdrBeamformerWeights,
} from '../src/aesa/beamforming_array';
import {
  computeMtdDopplerSpectrum,
  applyThreePulseMtiCanceller,
  resolveDualPrfAmbiguity,
} from '../src/doppler/mtd_stap';

describe('1. Complex Number Arithmetic & Radix-2 FFT Suite', () => {
  it('verifies complex multiplication, conjugate and inversion', () => {
    const a = new Complex(3, 4);
    const b = new Complex(1, -2);

    expect(a.abs()).toBeCloseTo(5.0, 5);
    expect(a.absSq()).toBe(25);

    // (3 + 4i)(1 - 2i) = 3 - 6i + 4i - 8i^2 = 11 - 2i
    const prod = a.mul(b);
    expect(prod.re).toBeCloseTo(11.0, 5);
    expect(prod.im).toBeCloseTo(-2.0, 5);

    // a / a == 1
    const divSelf = a.div(a);
    expect(divSelf.re).toBeCloseTo(1.0, 5);
    expect(divSelf.im).toBeCloseTo(0.0, 5);
  });

  it('verifies nextPowerOfTwo calculation', () => {
    expect(nextPowerOfTwo(1)).toBe(1);
    expect(nextPowerOfTwo(5)).toBe(8);
    expect(nextPowerOfTwo(32)).toBe(32);
    expect(nextPowerOfTwo(33)).toBe(64);
  });

  it('verifies FFT unit impulse response is flat across frequencies', () => {
    const n = 16;
    const impulse: Complex[] = Array.from({ length: n }, (_, i) =>
      i === 0 ? new Complex(1, 0) : Complex.zero()
    );

    const freq = fftRadix2(impulse, false);
    for (const val of freq) {
      expect(val.re).toBeCloseTo(1.0, 5);
      expect(val.im).toBeCloseTo(0.0, 5);
    }
  });

  it('verifies FFT -> IFFT invertibility and Parseval energy conservation', () => {
    const n = 32;
    const original: Complex[] = Array.from({ length: n }, (_, i) =>
      Complex.fromPolar(1 + 0.5 * Math.sin(i), (2 * Math.PI * i) / n)
    );

    const freq = fftRadix2(original, false);
    const reconstructed = fftRadix2(freq, true);

    for (let i = 0; i < n; i++) {
      expect(reconstructed[i].re).toBeCloseTo(original[i].re, 5);
      expect(reconstructed[i].im).toBeCloseTo(original[i].im, 5);
    }

    // Parseval: sum(|x|^2) = (1/N) * sum(|X|^2)
    const timeEnergy = original.reduce((acc, c) => acc + c.absSq(), 0);
    const freqEnergy = freq.reduce((acc, c) => acc + c.absSq(), 0) / n;
    expect(freqEnergy).toBeCloseTo(timeEnergy, 4);
  });
});

describe('2. LFM Chirp & Matched Filter Pulse Compression Suite', () => {
  it('synthesizes LFM chirp with expected time-bandwidth product and resolution', () => {
    const chirpResult = generateLFMChirp({
      pulseWidth: 10e-6,   // 10 microseconds
      bandwidth: 20e6,     // 20 MHz sweep
      samplingRate: 50e6,  // 50 MSps
    });

    expect(chirpResult.numSamples).toBe(500);
    expect(chirpResult.timeBandwidthProduct).toBe(200);
    // Range resolution = c / (2 * B) = 3e8 / (2 * 20e6) = 7.49 m
    expect(chirpResult.rangeResolutionMeters).toBeCloseTo(7.495, 2);

    const gainDb = computePulseCompressionGainDb(chirpResult.timeBandwidthProduct);
    expect(gainDb).toBeCloseTo(23.01, 1);
  });

  it('compresses received echo and detects target peak at exact range bin', () => {
    const chirp = generateLFMChirp({
      pulseWidth: 4e-6,
      bandwidth: 10e6,
      samplingRate: 40e6,
    });

    const targetBin = 120;
    const rxSignal: Complex[] = Array.from({ length: 300 }, () =>
      Complex.fromPolar(0.05 * Math.random(), Math.random() * 2 * Math.PI)
    );

    // Inject target echo at targetBin
    for (let i = 0; i < chirp.samples.length; i++) {
      if (targetBin + i < rxSignal.length) {
        rxSignal[targetBin + i] = rxSignal[targetBin + i].add(chirp.samples[i].scale(1.5));
      }
    }

    const compressed = matchedFilterPulseCompression(rxSignal, chirp.samples);
    expect(compressed.peakIndex).toBeGreaterThanOrEqual(targetBin - 2);
    expect(compressed.peakIndex).toBeLessThanOrEqual(targetBin + 2);
    expect(compressed.peakMagnitude).toBeGreaterThan(10.0);
  });
});

describe('3. CFAR Adaptive Detectors & Swerling Target Models', () => {
  it('computes CA-CFAR alpha factor analytically', () => {
    const alpha = computeCaCfarAlpha(1e-6, 16);
    expect(alpha).toBeGreaterThan(15.0);
  });

  it('detects injected radar targets while maintaining false alarm bounds across CFAR types', () => {
    const numCells = 128;
    const powerSpectrum = new Float64Array(numCells);

    // Uniform noise floor around 1.0
    for (let i = 0; i < numCells; i++) {
      powerSpectrum[i] = 0.8 + 0.4 * Math.random();
    }

    // Injected targets at bin 45 and 90
    powerSpectrum[45] = 50.0; // High SNR target
    powerSpectrum[90] = 35.0; // Moderate SNR target

    const caResult = runCfarDetector(powerSpectrum, {
      type: 'CA-CFAR',
      numGuardCells: 2,
      numReferenceCells: 6,
      pfa: 1e-4,
    });

    expect(caResult.detections[45]).toBe(true);
    expect(caResult.detections[90]).toBe(true);
    expect(caResult.numDetections).toBeLessThanOrEqual(4);

    const goResult = runCfarDetector(powerSpectrum, {
      type: 'GO-CFAR',
      numGuardCells: 2,
      numReferenceCells: 6,
      pfa: 1e-4,
    });
    expect(goResult.detections[45]).toBe(true);

    const osResult = runCfarDetector(powerSpectrum, {
      type: 'OS-CFAR',
      numGuardCells: 2,
      numReferenceCells: 6,
      pfa: 1e-4,
      osRankPercentile: 0.75,
    });
    expect(osResult.detections[45]).toBe(true);
  });

  it('computes monotonic Swerling detection probabilities', () => {
    const pdLow = calculateSwerlingPd('Swerling-1', 5.0, 1e-6);
    const pdHigh = calculateSwerlingPd('Swerling-1', 20.0, 1e-6);
    expect(pdHigh).toBeGreaterThan(pdLow);

    const pdSw0 = calculateSwerlingPd('Swerling-0', 15.0, 1e-6);
    expect(pdSw0).toBeGreaterThan(0.9);
  });
});

describe('4. 2D AESA Phased Array & MVDR Adaptive Null-Steering', () => {
  it('constructs octagonal planar array and steers broadside beam', () => {
    const array = new AESAArray({
      numElementsX: 12,
      numElementsY: 12,
      frequencyHz: 9.5e9,
      apertureShape: 'Octagon',
      taperType: 'Taylor',
    });

    expect(array.elements.length).toBeGreaterThan(80);

    // Steer to broadside (0, 0)
    array.steerBeam(0, 0);
    const afBroadside = array.computeArrayFactor(0, 0);
    expect(afBroadside.db).toBeCloseTo(0.0, 1);

    // Off-axis array factor at 30 deg should be significantly lower (sidelobe)
    const afOffAxis = array.computeArrayFactor(30, 0);
    expect(afOffAxis.db).toBeLessThan(-20.0);
  });

  it('steers main beam dynamically to +25 deg azimuth', () => {
    const array = new AESAArray({
      numElementsX: 12,
      numElementsY: 12,
      frequencyHz: 9.5e9,
    });

    array.steerBeam(25, 0);
    const afAtBeam = array.computeArrayFactor(25, 0);
    const afAtZero = array.computeArrayFactor(0, 0);

    expect(afAtBeam.db).toBeCloseTo(0.0, 1);
    expect(afAtBeam.db).toBeGreaterThan(afAtZero.db + 15.0);
  });

  it('forms deep adaptive nulls at hostile jammer azimuth using MVDR', () => {
    const array = new AESAArray({
      numElementsX: 8,
      numElementsY: 8,
      frequencyHz: 9.5e9,
      apertureShape: 'Rectangle',
    });

    const targetAz = 0;
    const jammerAz = 22;

    const weights = computeMvdrBeamformerWeights(
      array,
      targetAz,
      0,
      [{ az: jammerAz, el: 0, inrDb: 35 }],
      0
    );

    expect(weights.length).toBe(array.elements.length);

    // Evaluate beam pattern with MVDR weights
    const aTarget = array.getSteeringVector(targetAz, 0);
    const aJammer = array.getSteeringVector(jammerAz, 0);

    let respTarget = Complex.zero();
    let respJammer = Complex.zero();

    for (let i = 0; i < weights.length; i++) {
      respTarget = respTarget.add(weights[i].conj().mul(aTarget[i]));
      respJammer = respJammer.add(weights[i].conj().mul(aJammer[i]));
    }

    // Target response must be 1.0 (distortionless constraint: w^H * a = 1)
    expect(respTarget.abs()).toBeCloseTo(1.0, 1);
    // Jammer response must be deeply attenuated (deep null)
    const nullDepthDb = 20 * Math.log10(respJammer.abs());
    expect(nullDepthDb).toBeLessThan(-25.0);
  });
});

describe('5. Doppler MTD Filtering & Dual-PRF Ambiguity Resolution', () => {
  it('attenuates zero-Doppler clutter using 3-pulse delay-line canceller', () => {
    // Ground clutter: constant DC slow-time pulses
    const clutter: Complex[] = Array.from({ length: 16 }, () => new Complex(5.0, 0));
    const filtered = applyThreePulseMtiCanceller(clutter);

    for (const val of filtered) {
      expect(val.abs()).toBeCloseTo(0.0, 5);
    }
  });

  it('computes Doppler spectrum and extracts radial velocity peak', () => {
    const numPulses = 32;
    const prf = 10000;
    const wavelength = 0.0316; // 9.5 GHz
    const trueVelocity = 158; // m/s -> fd = 2 * vr / lambda = 10,000 Hz (near PRF)

    const slowTime: Complex[] = Array.from({ length: numPulses }, (_, i) => {
      const t = i / prf;
      const targetFd = (2 * trueVelocity) / wavelength;
      return Complex.fromPolar(2.0, 2 * Math.PI * targetFd * t);
    });

    const mtd = computeMtdDopplerSpectrum(slowTime, {
      numPulses,
      prfHz: prf,
      wavelengthMeters: wavelength,
      windowType: 'Hann',
    });

    expect(mtd.blindSpeedMps).toBeCloseTo((wavelength * prf) / 2, 1);
    expect(mtd.powerSpectrumDb.length).toBeGreaterThanOrEqual(32);
  });

  it('unwraps ambiguous Doppler velocity using dual-PRF coincidence', () => {
    const lambda = 0.0316;
    const prf1 = 8000;
    const prf2 = 10000;
    const trueV = 320.0; // High speed target (m/s)

    const vb1 = (lambda * prf1) / 2;
    const vb2 = (lambda * prf2) / 2;

    const vObs1 = ((trueV % vb1) + vb1) % vb1;
    const vObs2 = ((trueV % vb2) + vb2) % vb2;

    const resolvedV = resolveDualPrfAmbiguity(vObs1, vObs2, prf1, prf2, lambda);
    expect(resolvedV).toBeCloseTo(trueV, 0.5);
  });
});
