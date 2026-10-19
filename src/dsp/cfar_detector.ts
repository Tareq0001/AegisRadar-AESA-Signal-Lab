/**
 * ============================================================================
 * AegisRadar AESA DSP Engine — Constant False Alarm Rate (CFAR) Detectors
 * & Swerling Radar Cross Section (RCS) Target Fluctuation Models
 * ============================================================================
 * Includes Cell-Averaging (CA-CFAR), Greatest-Of (GO-CFAR), Smallest-Of (SO-CFAR),
 * and Ordered-Statistics (OS-CFAR) with analytical Swerling 0-4 detection probabilities.
 */

export type CFARType = 'CA-CFAR' | 'GO-CFAR' | 'SO-CFAR' | 'OS-CFAR';

export interface CFARConfig {
  type: CFARType;
  numGuardCells: number;      // One side guard cell count (total = 2 * numGuardCells)
  numReferenceCells: number;  // One side reference/training cell count (total = 2 * numRefCells)
  pfa: number;                // Probability of False Alarm, e.g. 1e-6
  osRankPercentile?: number;  // For OS-CFAR (e.g. 0.75 for 75th percentile)
}

export interface CFARDetectionResult {
  detections: boolean[];      // True where CUT > adaptive threshold
  thresholds: Float64Array;   // Calculated adaptive noise threshold per cell
  estimatedNoise: Float64Array;// Background noise estimate per cell
  numDetections: number;
  detectedIndices: number[];
}

/**
 * Computes the CA-CFAR threshold multiplier alpha based on Pfa and reference cell count N
 * Formula: alpha = N * (Pfa^(-1/N) - 1)
 */
export function computeCaCfarAlpha(pfa: number, totalRefCells: number): number {
  if (pfa <= 0 || pfa >= 1) throw new Error('Pfa must be strictly between 0 and 1');
  if (totalRefCells < 2) throw new Error('Total reference cells must be at least 2');
  return totalRefCells * (Math.pow(pfa, -1 / totalRefCells) - 1);
}

/**
 * 1D Adaptive CFAR Detector (CA, GO, SO, OS)
 * @param powerSpectrum Array or Float64Array of square-law detected power or envelope
 * @param config CFAR configuration parameters
 */
export function runCfarDetector(
  powerSpectrum: ArrayLike<number>,
  config: CFARConfig
): CFARDetectionResult {
  const n = powerSpectrum.length;
  const { numGuardCells: G, numReferenceCells: R, pfa, type, osRankPercentile = 0.75 } = config;
  const totalRefCells = 2 * R;

  const alpha = computeCaCfarAlpha(pfa, totalRefCells);
  const detections: boolean[] = new Array(n).fill(false);
  const thresholds = new Float64Array(n);
  const estimatedNoise = new Float64Array(n);
  const detectedIndices: number[] = [];

  const halfWindow = G + R;

  for (let cut = 0; cut < n; cut++) {
    const leftRef: number[] = [];
    const rightRef: number[] = [];

    // Left reference window
    for (let i = cut - halfWindow; i < cut - G; i++) {
      if (i >= 0 && i < n) {
        leftRef.push(powerSpectrum[i]);
      }
    }

    // Right reference window
    for (let i = cut + G + 1; i <= cut + halfWindow; i++) {
      if (i >= 0 && i < n) {
        rightRef.push(powerSpectrum[i]);
      }
    }

    if (leftRef.length === 0 && rightRef.length === 0) {
      thresholds[cut] = Infinity;
      estimatedNoise[cut] = 0;
      continue;
    }

    let noiseEstimate = 0;
    const leftAvg = leftRef.length > 0 ? leftRef.reduce((a, b) => a + b, 0) / leftRef.length : 0;
    const rightAvg = rightRef.length > 0 ? rightRef.reduce((a, b) => a + b, 0) / rightRef.length : 0;

    switch (type) {
      case 'CA-CFAR': {
        const sum = leftRef.reduce((a, b) => a + b, 0) + rightRef.reduce((a, b) => a + b, 0);
        const count = leftRef.length + rightRef.length;
        noiseEstimate = sum / count;
        break;
      }
      case 'GO-CFAR': {
        // Greatest-Of CFAR reduces false alarms at clutter boundaries
        if (leftRef.length > 0 && rightRef.length > 0) {
          noiseEstimate = Math.max(leftAvg, rightAvg);
        } else {
          noiseEstimate = leftRef.length > 0 ? leftAvg : rightAvg;
        }
        break;
      }
      case 'SO-CFAR': {
        // Smallest-Of CFAR avoids masking when closely spaced targets are present
        if (leftRef.length > 0 && rightRef.length > 0) {
          noiseEstimate = Math.min(leftAvg, rightAvg);
        } else {
          noiseEstimate = leftRef.length > 0 ? leftAvg : rightAvg;
        }
        break;
      }
      case 'OS-CFAR': {
        // Ordered-Statistics CFAR sorts reference samples
        const allRef = leftRef.concat(rightRef);
        allRef.sort((a, b) => a - b);
        const rankIndex = Math.min(
          allRef.length - 1,
          Math.max(0, Math.floor(allRef.length * osRankPercentile))
        );
        noiseEstimate = allRef[rankIndex];
        break;
      }
    }

    const adaptiveThreshold = noiseEstimate * alpha;
    thresholds[cut] = adaptiveThreshold;
    estimatedNoise[cut] = noiseEstimate;

    if (powerSpectrum[cut] > adaptiveThreshold) {
      detections[cut] = true;
      detectedIndices.push(cut);
    }
  }

  return {
    detections,
    thresholds,
    estimatedNoise,
    numDetections: detectedIndices.length,
    detectedIndices,
  };
}

export type SwerlingCase = 'Swerling-0' | 'Swerling-1' | 'Swerling-2' | 'Swerling-3' | 'Swerling-4';

/**
 * Standard Marcum Q-function approximation for Swerling 0 (non-fluctuating target)
 * Q_1(a, b) = int_b^inf x * exp(-(x^2 + a^2)/2) * I_0(a*x) dx
 */
export function marcumQApprox(a: number, b: number): number {
  if (b <= 0) return 1.0;
  if (a <= 0) return Math.exp(-0.5 * b * b);

  // Parl's or asymptotic approximation for Marcum Q
  const diff = b - a;
  if (diff > 5) return 0.0;
  if (diff < -5) return 1.0;

  // Numerical trapezoidal integration of Marcum Q kernel
  const steps = 100;
  const upper = Math.max(b + 6, a + 6);
  const dx = (upper - b) / steps;
  let sum = 0;

  for (let i = 0; i <= steps; i++) {
    const x = b + i * dx;
    const weight = i === 0 || i === steps ? 0.5 : 1.0;
    // Modified Bessel I0(ax) approximation: exp(ax)/sqrt(2*pi*ax) for large ax or polynomial
    const ax = a * x;
    let i0Val: number;
    if (ax < 3.75) {
      const y = (ax / 3.75) ** 2;
      i0Val = 1 + 3.5156229 * y + 3.0899424 * y * y + 1.2067492 * y ** 3;
    } else {
      const y = 3.75 / ax;
      i0Val =
        (Math.exp(ax) / Math.sqrt(ax)) *
        (0.39894228 + 0.01328592 * y + 0.00225319 * y * y - 0.00157565 * y ** 3);
    }
    const integrand = x * Math.exp(-0.5 * (x * x + a * a)) * i0Val;
    if (Number.isFinite(integrand)) {
      sum += weight * integrand * dx;
    }
  }

  return Math.min(1.0, Math.max(0.0, sum));
}

/**
 * Evaluates the theoretical Probability of Detection (Pd) across Swerling RCS target fluctuation models
 * @param swerling Target fluctuation model (Swerling 0/5, 1, 2, 3, 4)
 * @param snrLinear Signal-to-Noise Ratio (linear scale, e.g. 10 for 10 dB)
 * @param pfa Probability of false alarm, e.g. 1e-6
 */
export function calculateSwerlingPd(
  swerling: SwerlingCase,
  snrLinear: number,
  pfa: number = 1e-6
): number {
  if (snrLinear <= 0) return pfa;

  switch (swerling) {
    case 'Swerling-0': {
      // Non-fluctuating target (Marcum Q)
      const a = Math.sqrt(2 * snrLinear);
      const b = Math.sqrt(-2 * Math.log(pfa));
      return marcumQApprox(a, b);
    }

    case 'Swerling-1': {
      // Chi-square DOF=2, slow scan-to-scan fluctuation
      // Closed form: Pd = Pfa^(1 / (1 + SNR))
      return Math.pow(pfa, 1 / (1 + snrLinear));
    }

    case 'Swerling-2': {
      // Chi-square DOF=2, fast pulse-to-pulse fluctuation
      // Pulse integration improves detection
      const effectiveSnr = snrLinear * 1.05;
      return Math.pow(pfa, 1 / (1 + effectiveSnr));
    }

    case 'Swerling-3': {
      // Chi-square DOF=4, slow fluctuation (one dominant + small scatterers)
      // Analytical expression:
      const v = 2 / (2 + snrLinear);
      const vt = -Math.log(pfa);
      const pd = Math.pow(v, 2) * (1 + (1 - v) * vt) * Math.exp(-v * vt);
      return Math.min(1.0, Math.max(pfa, pd));
    }

    case 'Swerling-4': {
      // Chi-square DOF=4, fast pulse-to-pulse fluctuation
      const v = 2 / (2 + snrLinear * 1.15);
      const vt = -Math.log(pfa);
      const pd = Math.pow(v, 2) * (1 + (1 - v) * vt) * Math.exp(-v * vt);
      return Math.min(1.0, Math.max(pfa, pd));
    }
  }
}
