# AegisRadar AESA — 3D Phased Array Radar & Pulse Doppler DSP Lab

[![Live Demo](https://img.shields.io/badge/Live_AESA_Lab-Tareq0001.github.io-emerald?style=for-the-badge&logo=github)](https://tareq0001.github.io/AegisRadar-AESA-Signal-Lab/)
[![Frequency](https://img.shields.io/badge/Band-X--Band_9.5_GHz_%7C_GaN_T%2FR-blue?style=for-the-badge)](https://tareq0001.github.io/AegisRadar-AESA-Signal-Lab/)
[![Three.js](https://img.shields.io/badge/WebGL-Three.js_Beamforming_Engine-orange?style=for-the-badge&logo=three.js)](https://tareq0001.github.io/AegisRadar-AESA-Signal-Lab/)

> **A precision 3D Active Electronically Scanned Array (AESA) radar and digital signal processing (DSP) interactive laboratory. Features Gallium Nitride (GaN) transmit/receive module phase steering, digital beamforming in azimuth and elevation, an interactive real-time Radar Range Equation solver, Plan Position Indicator (PPI) polar radar scope, and Doppler FFT velocity spectrum with Moving Target Indication (MTI) ground clutter filtering.**

---

## ⚡ Core Radar & DSP Engineering Capabilities

1. **3D AESA Antenna Array & Digital Beamforming**
   - Octagonal aperture populated with hundreds of Gallium Nitride (GaN) Transmit/Receive (T/R) modules.
   - Independent phase shifting across elements simulating digital beam steering:
     $$\Delta \phi = \frac{2\pi}{\lambda} d \sin\theta$$
   - Visual 3D steered main beam radiation cone and dynamic phase color shading.

2. **Real-Time Radar Range Equation Mathematical Solver**
   - Computes received signal power ($P_r$) and Signal-to-Noise Ratio (SNR) continuously:
     $$P_r = \frac{P_t \cdot G^2 \cdot \lambda^2 \cdot \sigma}{(4\pi)^3 \cdot R^4}$$
   - Realistic parameters: $P_t = 120\text{ kW}$ peak power, $G = 42.5\text{ dBi}$ gain, $f = 9.5\text{ GHz}$ (X-Band wavelength $\lambda = 3.16\text{ cm}$).
   - Target Radar Cross Section (RCS - $\sigma$) presets: Stealth aircraft ($0.005\text{ m}^2$), Loitering Drone ($0.02\text{ m}^2$), Cruise missile ($0.1\text{ m}^2$), Fighter ($3.0\text{ m}^2$).

3. **Multi-Scope Radar Displays & Signal Processing**
   - **Plan Position Indicator (PPI Scope):** Polar radar sweep with range rings up to 180 km and real-time target blip detection.
   - **Doppler FFT Spectrum:** Fast Fourier Transform velocity spectrum calculating target radial velocity ($v_r = 612\text{ m/s}$) and Doppler shift ($f_d = +38.7\text{ kHz}$).
   - **Moving Target Indication (MTI):** High-pass clutter rejection filter isolating high-speed threats from background noise.

---

## 🛠️ Architecture & Tech Stack

- **Graphics:** Three.js (r128), Custom shader beam lobes, Canvas 2D PPI & FFT waterfall.
- **Electromagnetics:** AESA array phase calculations, Radar Range Equation, Doppler frequency shifts.
- **Audio:** Web Audio API procedural radar chirp & ping synthesis.
- **UI:** Tactile engineering deck, bilingual English/Arabic interface.

---

## 🚀 Live Deployment

- **Online AESA Radar Lab:** [https://tareq0001.github.io/AegisRadar-AESA-Signal-Lab/](https://tareq0001.github.io/AegisRadar-AESA-Signal-Lab/)
- **Author:** Eng. Tareq Aboushi (Tareq0001)
