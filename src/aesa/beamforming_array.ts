/**
 * ============================================================================
 * AegisRadar AESA Signal Lab — 2D Phased Array Digital Beamforming (DBF)
 * Taylor Tapering, 3D Array Factor & Adaptive MVDR Jammer Null-Steering
 * ============================================================================
 * Models enterprise Active Electronically Scanned Arrays (AESA) with independent
 * GaN Transmit/Receive (T/R) phase and amplitude control.
 */

import { Complex } from '../dsp/fft_chirp';

export interface ArrayElement {
  id: number;
  x: number;          // Position along X axis in meters
  y: number;          // Position along Y axis in meters
  z: number;          // Position along Z axis in meters
  amplitude: number;  // Tapering amplitude weight w_n
  phaseRad: number;   // Steering phase shift in radians
}

export type TaperWindowType = 'Uniform' | 'Taylor' | 'Hamming' | 'Chebyshev';

export interface AESAConfig {
  numElementsX: number;     // e.g. 16
  numElementsY: number;     // e.g. 16
  frequencyHz: number;      // e.g. 9.5e9 (9.5 GHz X-Band)
  elementSpacingRatio?: number; // d / lambda, default 0.5 (half-wavelength)
  taperType?: TaperWindowType;
  taylorSLLDb?: number;     // e.g. -30 dB sidelobe level
  taylorNBar?: number;      // e.g. 4
  apertureShape?: 'Rectangle' | 'Octagon' | 'Circle';
}

const SPEED_OF_LIGHT = 299792458; // m/s

/**
 * 2D Planar AESA Phased Array Model
 */
export class AESAArray {
  public readonly elements: ArrayElement[] = [];
  public readonly wavelength: number;
  public readonly k: number; // Wave number 2*pi / lambda
  public readonly dx: number;
  public readonly dy: number;
  public steeredAzDeg: number = 0;
  public steeredElDeg: number = 0;

  constructor(public readonly config: AESAConfig) {
    this.wavelength = SPEED_OF_LIGHT / config.frequencyHz;
    this.k = (2 * Math.PI) / this.wavelength;
    const spacingRatio = config.elementSpacingRatio ?? 0.5;
    this.dx = spacingRatio * this.wavelength;
    this.dy = spacingRatio * this.wavelength;

    this.buildArrayElements();
    this.applyAmplitudeTaper(config.taperType ?? 'Taylor');
    this.steerBeam(0, 0);
  }

  private buildArrayElements(): void {
    const { numElementsX: nx, numElementsY: ny, apertureShape = 'Octagon' } = this.config;
    const x0 = ((nx - 1) * this.dx) / 2;
    const y0 = ((ny - 1) * this.dy) / 2;
    const maxR = Math.max(x0, y0);

    let id = 0;
    for (let iy = 0; iy < ny; iy++) {
      for (let ix = 0; ix < nx; ix++) {
        const x = ix * this.dx - x0;
        const y = iy * this.dy - y0;

        // Apply aperture boundary filter
        if (apertureShape === 'Circle') {
          if (Math.hypot(x, y) > maxR * 1.05) continue;
        } else if (apertureShape === 'Octagon') {
          // Octagonal chamfer: |x| + |y| <= 1.4 * maxR
          if (Math.abs(x) + Math.abs(y) > 1.35 * maxR) continue;
        }

        this.elements.push({
          id: id++,
          x,
          y,
          z: 0,
          amplitude: 1.0,
          phaseRad: 0.0,
        });
      }
    }
  }

  /**
   * Applies amplitude tapering across the array to suppress sidelobes
   */
  public applyAmplitudeTaper(taper: TaperWindowType): void {
    const maxDist = Math.max(...this.elements.map(e => Math.hypot(e.x, e.y)));

    for (const el of this.elements) {
      const rNorm = maxDist > 0 ? Math.hypot(el.x, el.y) / maxDist : 0;

      switch (taper) {
        case 'Uniform':
          el.amplitude = 1.0;
          break;
        case 'Hamming':
          el.amplitude = 0.54 + 0.46 * Math.cos(Math.PI * rNorm);
          break;
        case 'Taylor':
        case 'Chebyshev': {
          // Circular Taylor n-bar approximation taper
          // Tapers from 1.0 at center down to ~0.25 at aperture edge for ~30 dB SLL
          const edgePedestal = 0.28;
          el.amplitude = edgePedestal + (1.0 - edgePedestal) * Math.pow(1 - Math.min(1, rNorm * rNorm), 1.2);
          break;
        }
      }
    }
  }

  /**
   * Steers the main beam to given Azimuth and Elevation (in degrees)
   */
  public steerBeam(azimuthDeg: number, elevationDeg: number): void {
    this.steeredAzDeg = azimuthDeg;
    this.steeredElDeg = elevationDeg;

    const azRad = (azimuthDeg * Math.PI) / 180;
    const elRad = (elevationDeg * Math.PI) / 180;

    // Direction cosines: u = sin(az)*cos(el), v = sin(el)
    const u0 = Math.sin(azRad) * Math.cos(elRad);
    const v0 = Math.sin(elRad);

    for (const el of this.elements) {
      // Phase delay: phi_n = -k * (x_n * u0 + y_n * v0)
      el.phaseRad = -this.k * (el.x * u0 + el.y * v0);
    }
  }

  /**
   * Computes the normalized Array Factor (AF) in dB at observation angle (theta, phi)
   */
  public computeArrayFactor(azimuthDeg: number, elevationDeg: number): {
    linear: number;
    db: number;
    phase: number;
  } {
    const azRad = (azimuthDeg * Math.PI) / 180;
    const elRad = (elevationDeg * Math.PI) / 180;

    const u = Math.sin(azRad) * Math.cos(elRad);
    const v = Math.sin(elRad);

    let sumRe = 0;
    let sumIm = 0;
    let maxPossibleAmp = 0;

    for (const el of this.elements) {
      maxPossibleAmp += el.amplitude;
      // Net phase: phase steering + spatial propagation phase
      const phase = el.phaseRad + this.k * (el.x * u + el.y * v);
      sumRe += el.amplitude * Math.cos(phase);
      sumIm += el.amplitude * Math.sin(phase);
    }

    const linear = Math.hypot(sumRe, sumIm) / (maxPossibleAmp || 1);
    const db = 20 * Math.log10(Math.max(1e-6, linear));
    const phase = Math.atan2(sumIm, sumRe);

    return { linear, db, phase };
  }

  /**
   * Generates a spatial steering vector a(az, el) for the array
   */
  public getSteeringVector(azimuthDeg: number, elevationDeg: number): Complex[] {
    const azRad = (azimuthDeg * Math.PI) / 180;
    const elRad = (elevationDeg * Math.PI) / 180;
    const u = Math.sin(azRad) * Math.cos(elRad);
    const v = Math.sin(elRad);

    return this.elements.map(el => {
      const spatialPhase = this.k * (el.x * u + el.y * v);
      return Complex.fromPolar(1.0, spatialPhase);
    });
  }
}

/**
 * Minimum Variance Distortionless Response (MVDR) / Capon Adaptive Beamformer
 * Synthesizes adaptive element weights to form deep nulls in jammer directions
 * while preserving unit response along target direction:
 * w_MVDR = (R_x^-1 * a(theta_0)) / (a(theta_0)^H * R_x^-1 * a(theta_0))
 */
export function computeMvdrBeamformerWeights(
  array: AESAArray,
  targetAz: number,
  targetEl: number,
  jammers: Array<{ az: number; el: number; inrDb: number }>,
  noisePowerDb: number = 0
): Complex[] {
  const N = array.elements.length;
  const aTarget = array.getSteeringVector(targetAz, targetEl);

  // Initialize Covariance Matrix R (N x N)
  const noisePower = Math.pow(10, noisePowerDb / 10);
  const R: Complex[][] = Array.from({ length: N }, (_, i) =>
    Array.from({ length: N }, (_, j) => (i === j ? new Complex(noisePower, 0) : Complex.zero()))
  );

  // Add interference sources to covariance matrix: R += sum(P_j * a_j * a_j^H)
  for (const jammer of jammers) {
    const jammerPower = Math.pow(10, jammer.inrDb / 10);
    const aJammer = array.getSteeringVector(jammer.az, jammer.el);

    for (let i = 0; i < N; i++) {
      for (let j = 0; j < N; j++) {
        // a_j[i] * conj(a_j[j])
        const term = aJammer[i].mul(aJammer[j].conj()).scale(jammerPower);
        R[i][j] = R[i][j].add(term);
      }
    }
  }

  // Regularize diagonal loading for matrix inversion numerical stability
  const diagLoad = noisePower * 1e-4;
  for (let i = 0; i < N; i++) {
    R[i][i] = R[i][i].add(new Complex(diagLoad, 0));
  }

  // Invert Covariance Matrix R using Gauss-Jordan elimination
  const RInv = invertComplexMatrix(R);

  // Compute RInv * aTarget
  const rInvATarget: Complex[] = new Array(N);
  for (let i = 0; i < N; i++) {
    let sum = Complex.zero();
    for (let j = 0; j < N; j++) {
      sum = sum.add(RInv[i][j].mul(aTarget[j]));
    }
    rInvATarget[i] = sum;
  }

  // Denominator: aTarget^H * RInv * aTarget
  let denom = Complex.zero();
  for (let i = 0; i < N; i++) {
    denom = denom.add(aTarget[i].conj().mul(rInvATarget[i]));
  }

  // w_MVDR = rInvATarget / denom
  return rInvATarget.map(val => val.div(denom));
}

/**
 * Performs Gauss-Jordan elimination to invert an N x N complex matrix
 */
function invertComplexMatrix(mat: Complex[][]): Complex[][] {
  const n = mat.length;
  // Augment with identity matrix
  const A: Complex[][] = mat.map((row, i) => [
    ...row.map(c => new Complex(c.re, c.im)),
    ...Array.from({ length: n }, (_, j) => (i === j ? new Complex(1, 0) : Complex.zero())),
  ]);

  for (let col = 0; col < n; col++) {
    // Find pivot
    let maxRow = col;
    let maxMag = A[col][col].abs();
    for (let r = col + 1; r < n; r++) {
      const mag = A[r][col].abs();
      if (mag > maxMag) {
        maxMag = mag;
        maxRow = r;
      }
    }

    if (maxMag < 1e-12) {
      throw new Error('Singular matrix in MVDR inversion');
    }

    // Swap rows
    if (maxRow !== col) {
      const temp = A[col];
      A[col] = A[maxRow];
      A[maxRow] = temp;
    }

    // Scale pivot row
    const pivot = A[col][col];
    for (let c = col; c < 2 * n; c++) {
      A[col][c] = A[col][c].div(pivot);
    }

    // Eliminate other rows
    for (let r = 0; r < n; r++) {
      if (r !== col) {
        const factor = A[r][col];
        for (let c = col; c < 2 * n; c++) {
          A[r][c] = A[r][c].sub(factor.mul(A[col][c]));
        }
      }
    }
  }

  // Extract right half (inverse)
  return A.map(row => row.slice(n));
}
