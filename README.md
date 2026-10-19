# AegisRadar AESA — 3D Phased Array Radar & Coherent DSP Lab

[![CI Pipeline](https://github.com/Tareq0001/AegisRadar-AESA-Signal-Lab/actions/workflows/ci.yml/badge.svg)](https://github.com/Tareq0001/AegisRadar-AESA-Signal-Lab/actions/workflows/ci.yml)
[![Live Demo](https://img.shields.io/badge/Live_AESA_Lab-tareq0001.github.io-emerald?style=for-the-badge&logo=github)](https://tareq0001.github.io/AegisRadar-AESA-Signal-Lab/)
[![Band](https://img.shields.io/badge/Radar_Band-X--Band_9.5_GHz_%7C_GaN_T%2FR-blue?style=for-the-badge)](https://tareq0001.github.io/AegisRadar-AESA-Signal-Lab/)
[![Stack](https://img.shields.io/badge/Languages-TypeScript_%7C_C%2B%2B17_%7C_Python-teal?style=for-the-badge)](https://github.com/Tareq0001/AegisRadar-AESA-Signal-Lab)
[![Three.js](https://img.shields.io/badge/WebGL-Three.js_3D_Engine-orange?style=for-the-badge&logo=three.js)](https://threejs.org/)
[![License](https://img.shields.io/badge/License-MIT-green?style=for-the-badge)](LICENSE)

> **Enterprise-grade 3D Active Electronically Scanned Array (AESA) radar simulation, digital signal processing (DSP), adaptive beamforming, and electronic protection laboratory. Implements real-time Gallium Nitride (GaN) phase steering, Taylor amplitude tapering, Minimum Variance Distortionless Response (MVDR) jammer null-steering, Linear Frequency Modulation (LFM) matched filtering, Constant False Alarm Rate (CFAR) detection across Swerling RCS models, Moving Target Detection (MTD), and Space-Time Adaptive Processing (STAP).**

---

## 🛰️ Interactive WebGL Laboratory

Access the fully operational radar engineering console directly in your browser:  
🔗 **[https://tareq0001.github.io/AegisRadar-AESA-Signal-Lab/](https://tareq0001.github.io/AegisRadar-AESA-Signal-Lab/)**

- **3D AESA Physical Array Aperture:** Octagonal planar array with hundreds of individual GaN T/R modules with real-time phase shading and radiation cone orientation.
- **Dynamic 2D Radar Scopes:** Plan Position Indicator (PPI) polar radar scope (180 km instrumented range) and Doppler Fast Fourier Transform (FFT) velocity waterfall.
- **Tactile Electronic Warfare Deck:** Azimuth ($\pm 60^\circ$) and Elevation ($-20^\circ \text{ to } +60^\circ$) beam steering sliders, RCS threat presets, and procedural RF pulse audio synthesis.
- **Bilingual Interface:** Full English & Arabic localization.

---

## 🏛️ System Architecture & Signal Chain

```
                +-------------------------------------------------+
                |      AESA 2D Planar Array Aperture (GaN T/R)    |
                +-------------------------------------------------+
                                        |
                 [Digital Phase Shifting: Δφ = -k(x·u + y·v)]
                 [Taylor Amplitude Windowing / Chebyshev Tapering]
                                        v
                +-------------------------------------------------+
                |        Digital Beamforming (DBF) & MVDR         |
                |   - Steers Mainlobe to (θ_az, θ_el)             |
                |   - Adaptive Null-Steering at Hostile Jammers   |
                +-------------------------------------------------+
                                        |
                                        v
                +-------------------------------------------------+
                |   Pulse Compression & Matched Filter Engine     |
                |   - Analytic LFM Chirp: s(t) = exp(j·π·K·t²)    |
                |   - Fast Frequency Convolution: Y(f) = X(f)·H*(f)|
                |   - Compression Gain: PCR = B · τ (~23.0 dB)    |
                +-------------------------------------------------+
                                        |
                                        v
                +-------------------------------------------------+
                |     Moving Target Detector (MTD) & STAP         |
                |   - 3-Pulse Clutter Delay-Line Canceller (MTI)  |
                |   - Doppler FFT Filter Bank across Pulses (CPI) |
                |   - Space-Time Clutter Ridge Suppression (STAP) |
                +-------------------------------------------------+
                                        |
                                        v
                +-------------------------------------------------+
                |          Adaptive CFAR Detection Suite          |
                |   - CA-CFAR: Homogeneous Thermal Background     |
                |   - GO-CFAR: Clutter Edge Discontinuities       |
                |   - SO-CFAR: Closely Spaced Target Resolution   |
                |   - OS-CFAR: Ordered-Statistics Ranking         |
                |   - Swerling 0-4 Fluctuation Probability Models |
                +-------------------------------------------------+
                                        |
                                        v
                +-------------------------------------------------+
                |    Track Extraction & Monopulse Estimation      |
                |   - Sum (Σ) & Difference (Δ_az, Δ_el) Channels  |
                |   - Tactical PPI Scope & Doppler Telemetry Deck |
                +-------------------------------------------------+
```

---

## 📐 Mathematical Formulations & DSP Theory

### 1. 2D AESA Planar Array Steering & Digital Beamforming
For an array of radiating elements located at coordinates $(x_n, y_n)$, the steering phase delay required to steer the electromagnetic mainlobe to spherical coordinates (Azimuth $\theta_{az}$, Elevation $\theta_{el}$) is:
$$u = \sin\theta_{az} \cos\theta_{el}, \quad v = \sin\theta_{el}, \quad k = \frac{2\pi}{\lambda}$$
$$\phi_n = -k (x_n u + y_n v)$$

The 3D normalized Array Factor (AF) including element amplitude weights $w_n$ is:
$$\text{AF}(\theta, \phi) = \frac{1}{\sum |w_n|} \left| \sum_{n=1}^{N} w_n e^{j (\phi_n + k(x_n \sin\theta \cos\phi + y_n \sin\theta \sin\phi))} \right|$$

### 2. Adaptive MVDR (Capon) Jammer Null-Steering
To neutralize hostile electronic countermeasure (ECM) stand-off jammers while preserving 0 dB gain along the target direction:
$$\vec{w}_{\text{MVDR}} = \frac{\mathbf{R}_x^{-1} \vec{a}(\theta_{\text{target}})}{\vec{a}^H(\theta_{\text{target}}) \mathbf{R}_x^{-1} \vec{a}(\theta_{\text{target}})}$$
Where $\mathbf{R}_x = \sigma_n^2 \mathbf{I} + \sum_{j} P_j \vec{a}(\theta_j) \vec{a}^H(\theta_j)$ is the sample covariance matrix. This synthesizes adaptive spatial nulls exceeding $-40\text{ dB}$ in the direction of jammers.

### 3. Radar Range Equation & Link Budget
Received signal echo power $P_r$ and signal-to-noise ratio (SNR) for a target with radar cross-section $\sigma$ at range $R$:
$$P_r = \frac{P_t \cdot G^2 \cdot \lambda^2 \cdot \sigma}{(4\pi)^3 \cdot R^4}$$
$$\text{SNR} = \frac{P_r}{k \cdot T_0 \cdot B \cdot F}$$
- Peak Transmit Power ($P_t$): $120\text{ kW}$
- Antenna Gain ($G$): $42.5\text{ dBi}$ ($\approx 17,782\text{ linear}$)
- Frequency ($f_0$): $9.5\text{ GHz}$ ($\lambda \approx 3.16\text{ cm}$)
- Receiver Noise Figure ($F$): $3.2\text{ dB}$
- Target RCS Presets ($\sigma$):
  - **Stealth Fighter:** $0.005\text{ m}^2$ (VLO signature)
  - **Loitering Drone / UAS:** $0.02\text{ m}^2$
  - **Supersonic Cruise Missile:** $0.1\text{ m}^2$
  - **Non-Stealth 4th-Gen Fighter:** $3.0\text{ m}^2$

### 4. LFM Chirp & Matched Filter Pulse Compression
An analytic linear frequency modulated (LFM) chirp pulse of duration $\tau$ and bandwidth $B$:
$$s(t) = \text{rect}\left(\frac{t}{\tau}\right) \exp\left(j \pi \frac{B}{\tau} t^2\right)$$
- **Pulse Compression Ratio (PCR):** $\text{PCR} = B \cdot \tau = 20\text{ MHz} \times 10\ \mu\text{s} = 200$
- **Theoretical Processing Gain:** $10 \log_{10}(B \cdot \tau) \approx 23.01\text{ dB}$
- **Compressed Range Resolution:** $\Delta R = \frac{c}{2B} = \frac{3 \times 10^8}{2 \times 20 \times 10^6} \approx 7.5\text{ m}$ (vs. $1.5\text{ km}$ uncompressed).

### 5. Constant False Alarm Rate (CFAR) Detectors
In non-stationary clutter, adaptive thresholding maintains constant false alarm probability $P_{fa}$:
- **CA-CFAR Multiplier:** $\alpha = N \left( P_{fa}^{-1/N} - 1 \right)$
- **Threshold:** $T = \alpha \cdot \hat{P}_n$
- **GO-CFAR:** $\hat{P}_n = \max(\bar{P}_{\text{left}}, \bar{P}_{\text{right}})$ (suppresses false alarms at clutter boundaries)
- **SO-CFAR:** $\hat{P}_n = \min(\bar{P}_{\text{left}}, \bar{P}_{\text{right}})$ (resolves closely spaced targets)
- **OS-CFAR:** $\hat{P}_n = P_{(k)}$ where $P_{(k)}$ is the $k$-th ordered statistic (robust against multiple interferers).

### 6. Doppler Processing, MTD & STAP
Target radial velocity $v_r$ produces a coherent Doppler frequency shift:
$$f_d = \frac{2 v_r}{\lambda}, \quad v_{\text{blind}} = \frac{\lambda \cdot \text{PRF}}{2}$$
- **3-Pulse Delay Line Canceller (MTI):** $H(z) = (1 - z^{-1})^2 = 1 - 2z^{-1} + z^{-2}$ (eliminates zero-Doppler stationary clutter).
- **Space-Time Adaptive Processing (STAP):** Suppresses airborne platform ground clutter ridges in 2D angle-Doppler space using 2D space-time steering vectors $\vec{v}(\theta, f_d) = \vec{a}(\theta) \otimes \vec{b}(f_d)$.

---

## 🗂️ Repository Structure

```
AegisRadar-AESA-Signal-Lab/
├── .github/
│   └── workflows/
│       └── ci.yml                     # Multi-language CI (TypeScript + C++17 + Python)
├── cpp_dsp/
│   ├── aesa_pulse_doppler.hpp         # C++17 SIMD Pulse-Doppler & Monopulse header
│   └── aesa_pulse_doppler.cpp         # C++17 Radix-2 FFT, matched filter & CA-CFAR
├── python_sim/
│   ├── radar_waveform_analyzer.py     # Woodward Ambiguity, Radar Range Eq & Swerling ROC
│   └── requirements.txt               # numpy, scipy, matplotlib
├── src/
│   ├── aesa/
│   │   └── beamforming_array.ts       # 2D Phased Array, Taylor taper & MVDR nulling
│   ├── doppler/
│   │   └── mtd_stap.ts                # Moving Target Detector & STAP clutter rejection
│   ├── dsp/
│   │   ├── cfar_detector.ts           # CA, GO, SO, OS-CFAR & Swerling 0-4 models
│   │   └── fft_chirp.ts               # Radix-2 FFT, LFM chirp & matched pulse compression
│   └── main.ts                        # Three.js 3D coordinator, Web Audio & radar UI
├── tests/
│   └── radar_dsp.test.ts              # Vitest DSP algorithmic verification suite
├── index.html                         # Modular HTML application entry point
├── package.json                       # Scripts, dependencies & metadata
├── tsconfig.json                      # Strict ES2022 TypeScript configuration
├── vite.config.ts                     # Bundler and development server configuration
└── README.md                          # Comprehensive technical documentation
```

---

## 🛠️ Build, Test & Execution Guide

### 1. TypeScript & Web Application
```bash
# Install dependencies
npm install

# Run static typechecking and Vitest algorithmic suite
npm run lint
npm run test

# Start local interactive development server
npm run dev
# -> Opens http://localhost:3000

# Build production bundle
npm run build
```

### 2. C++17 Pulse-Doppler & Monopulse Kernel
Compile the high-performance C++17 signal processing unit with maximum optimization:
```bash
# Linux / macOS / WSL
g++ -std=c++17 -O3 -Wall -Wextra -pedantic -c cpp_dsp/aesa_pulse_doppler.cpp -o cpp_dsp/aesa_pulse_doppler.o

# Or using package.json script
npm run build:cpp
```

### 3. Python Ambiguity Function & Waveform Analyzer
Execute the scientific waveform and radar range equation verification suite:
```bash
# Install scientific packages
pip install -r python_sim/requirements.txt

# Run automated algorithmic verification suite (exits 0 on 100% pass)
python python_sim/radar_waveform_analyzer.py --self-test
```

Sample output:
```text
============================================================================
  AEGISRADAR AESA SIGNAL LAB — RADAR WAVEFORM VERIFICATION SUITE
============================================================================
 [PASS] LFM Chirp Synthesized: 500 samples @ 50 MSps (Tau=10 µs, B=20 MHz)
 [PASS] Matched Filter Pulse Compression: TB=200, Gain=23.01 dB, Delta_R=8.99 m
 [PASS] Radar Range Equation: 50 km -> SNR=-4.4 dB, 100 km -> SNR=-16.4 dB
 [PASS] Swerling Detection Models: Pd(Sw0@15dB)=0.996, Pd(Sw1@15dB)=0.655
 [PASS] Woodward Ambiguity Surface: Peak=1.0000 centered at tau=0, nu=0
============================================================================
  ALL AEGISRADAR WAVEFORM & DSP ALGORITHMIC VERIFICATIONS PASSED (100%)
============================================================================
```

---

## 🌐 Production Deployment

- **Hosting:** GitHub Pages
- **Branch:** `gh-pages` (serves the pre-compiled standalone production bundle)
- **Live URL:** [https://tareq0001.github.io/AegisRadar-AESA-Signal-Lab/](https://tareq0001.github.io/AegisRadar-AESA-Signal-Lab/)

---

## 👨‍💻 Author & Attribution

**Eng. Tareq Aboushi** ([Tareq0001](https://github.com/Tareq0001))  
*AESA Radar Signal Processing, Digital Beamforming & Electronic Warfare Systems*

Licensed under the [MIT License](LICENSE).
