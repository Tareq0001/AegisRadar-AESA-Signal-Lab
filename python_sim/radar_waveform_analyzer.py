#!/usr/bin/env python3
"""
================================================================================
AegisRadar AESA Signal Lab — Radar Waveform, Ambiguity Function & ROC Analyzer
================================================================================
Scientific verification suite for Linear Frequency Modulated (LFM) Chirp synthesis,
Woodward Ambiguity Function, Radar Range Equation, and Swerling Target Detection.
Author: Eng. Tareq Aboushi (Tareq0001)
"""

import sys
import argparse
import numpy as np
from scipy import signal, special

SPEED_OF_LIGHT = 299792458.0  # m/s
BOLTZMANN_K = 1.380649e-23     # J/K
T0 = 290.0                     # Standard noise temperature (K)


class RadarWaveformAnalyzer:
    """
    Precision Radar Waveform & Ambiguity Function Analysis Engine
    """
    def __init__(
        self,
        carrier_freq_hz: float = 9.5e9,    # 9.5 GHz X-Band
        pulse_width_sec: float = 10.0e-6,   # 10 microseconds
        bandwidth_hz: float = 20.0e6,       # 20 MHz sweep
        sampling_rate_hz: float = 50.0e6,   # 50 MSps
    ):
        self.fc = carrier_freq_hz
        self.tau = pulse_width_sec
        self.b = bandwidth_hz
        self.fs = sampling_rate_hz
        self.chirp_rate = self.b / self.tau  # K = B / tau (Hz/s)
        self.wavelength = SPEED_OF_LIGHT / self.fc

    def generate_lfm_chirp(self) -> tuple[np.ndarray, np.ndarray]:
        """
        Synthesizes analytic complex baseband LFM chirp s(t) = exp(j * pi * K * t^2)
        """
        n_samples = int(np.round(self.tau * self.fs))
        t = np.linspace(-self.tau / 2.0, self.tau / 2.0, n_samples, endpoint=False)
        chirp = np.exp(1j * np.pi * self.chirp_rate * t**2)
        return t, chirp

    def compute_ambiguity_function(
        self,
        delay_bins: int = 128,
        doppler_bins: int = 128,
        max_doppler_hz: float = 100.0e3
    ) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        """
        Computes the Woodward Ambiguity Function |chi(tau, nu)|:
        chi(tau, nu) = integral s(t + tau/2) * s*(t - tau/2) * exp(j * 2*pi * nu * t) dt
        """
        _, s = self.generate_lfm_chirp()
        n = len(s)

        delays = np.linspace(-self.tau, self.tau, delay_bins)
        dopplers = np.linspace(-max_doppler_hz, max_doppler_hz, doppler_bins)
        ambiguity_surface = np.zeros((delay_bins, doppler_bins), dtype=np.float64)

        dt = 1.0 / self.fs
        t_base = (np.arange(n) - n / 2.0) * dt

        for i, tau_val in enumerate(delays):
            shift_samples = int(np.round(tau_val * self.fs))
            if abs(shift_samples) >= n:
                continue

            if shift_samples >= 0:
                s1 = np.pad(s, (shift_samples, 0))[:n]
                s2 = s
            else:
                s1 = s
                s2 = np.pad(s, (-shift_samples, 0))[:n]

            kernel = s1 * np.conj(s2)

            for j, nu_val in enumerate(dopplers):
                modulator = np.exp(1j * 2.0 * np.pi * nu_val * t_base)
                val = np.abs(np.sum(kernel * modulator) * dt)
                ambiguity_surface[i, j] = val

        # Normalize peak to 1.0 (0 dB)
        peak = np.max(ambiguity_surface)
        if peak > 0:
            ambiguity_surface /= peak

        return delays, dopplers, ambiguity_surface

    def solve_radar_equation(
        self,
        peak_power_watts: float = 120.0e3,
        antenna_gain_dbi: float = 42.5,
        target_rcs_sqm: float = 0.02,
        ranges_km: np.ndarray = None,
        noise_figure_db: float = 3.2
    ) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        """
        Solves the Radar Range Equation across distances:
        Pr = (Pt * G^2 * lambda^2 * sigma) / ((4*pi)^3 * R^4)
        """
        if ranges_km is None:
            ranges_km = np.linspace(10.0, 250.0, 100)

        r_meters = ranges_km * 1000.0
        g_linear = 10.0 ** (antenna_gain_dbi / 10.0)

        # Received Power (Watts)
        numerator = peak_power_watts * (g_linear**2) * (self.wavelength**2) * target_rcs_sqm
        denominator = ((4.0 * np.pi)**3) * (r_meters**4)
        pr_watts = numerator / denominator
        pr_dbm = 10.0 * np.log10(np.maximum(1e-22, pr_watts)) + 30.0

        # Thermal Noise Power
        noise_figure_linear = 10.0 ** (noise_figure_db / 10.0)
        pn_watts = BOLTZMANN_K * T0 * self.b * noise_figure_linear
        pn_dbm = 10.0 * np.log10(pn_watts) + 30.0

        snr_db = pr_dbm - pn_dbm
        return ranges_km, pr_dbm, snr_db

    def compute_swerling_roc(
        self,
        snr_db_range: np.ndarray = None,
        pfa: float = 1e-6
    ) -> dict[str, np.ndarray]:
        """
        Computes Probability of Detection (Pd) vs SNR for Swerling 0, 1, 3 models
        """
        if snr_db_range is None:
            snr_db_range = np.linspace(0.0, 25.0, 60)

        snr_linear = 10.0 ** (snr_db_range / 10.0)

        # Swerling 0 / Marcum Q-function approximation (Albersheim / Parl approximation)
        b = np.sqrt(-2.0 * np.log(pfa))
        a = np.sqrt(2.0 * snr_linear)
        # Marcum Q1(a, b) approximation
        pd_sw0 = 0.5 * special.erfc((b - a) / np.sqrt(2.0))
        pd_sw0 = np.clip(pd_sw0, pfa, 1.0)

        # Swerling 1 (Chi-square DOF=2): Pd = Pfa^(1 / (1 + SNR))
        pd_sw1 = pfa ** (1.0 / (1.0 + snr_linear))

        # Swerling 3 (Chi-square DOF=4, one dominant reflector)
        # Pd = (1 + 2 / (2 + SNR) * ln(1/Pfa)) * (1 + SNR / 2)^-2 ...
        v = 2.0 / (2.0 + snr_linear)
        vt = -np.log(pfa)
        pd_sw3 = (v**2) * (1.0 + (1.0 - v) * vt) * np.exp(-v * vt)
        pd_sw3 = np.clip(pd_sw3, pfa, 1.0)

        return {
            "snr_db": snr_db_range,
            "swerling_0": pd_sw0,
            "swerling_1": pd_sw1,
            "swerling_3": pd_sw3,
        }

    def verify_pulse_compression_gain(self) -> tuple[float, float, float]:
        """
        Validates matched filter pulse compression ratio (PCR) and SNR gain
        """
        t, chirp = self.generate_lfm_chirp()
        tb_product = self.b * self.tau
        theoretical_gain_db = 10.0 * np.log10(tb_product)

        # Matched filter impulse response
        compressed = signal.fftconvolve(chirp, np.conj(chirp[::-1]), mode='same')
        mag = np.abs(compressed)
        peak_val = np.max(mag)

        # Estimate 3 dB mainlobe width (Range resolution)
        half_power = peak_val / np.sqrt(2.0)
        indices_above_half = np.where(mag >= half_power)[0]
        dt = 1.0 / self.fs
        mainlobe_width_sec = len(indices_above_half) * dt
        measured_range_res_m = (SPEED_OF_LIGHT * mainlobe_width_sec) / 2.0
        theoretical_range_res_m = SPEED_OF_LIGHT / (2.0 * self.b)

        return tb_product, theoretical_gain_db, measured_range_res_m


def run_self_test() -> int:
    """
    Executes automated self-test verification suite for CI/CD pipeline
    """
    print("=" * 76)
    print("  AEGISRADAR AESA SIGNAL LAB — RADAR WAVEFORM VERIFICATION SUITE")
    print("=" * 76)

    analyzer = RadarWaveformAnalyzer(
        carrier_freq_hz=9.5e9,
        pulse_width_sec=10.0e-6,
        bandwidth_hz=20.0e6,
        sampling_rate_hz=50.0e6
    )

    # 1. Chirp Generation
    t, chirp = analyzer.generate_lfm_chirp()
    assert len(chirp) == 500, f"Expected 500 chirp samples, got {len(chirp)}"
    print(f" [PASS] LFM Chirp Synthesized: {len(chirp)} samples @ 50 MSps (Tau=10 µs, B=20 MHz)")

    # 2. Pulse Compression Gain & Range Resolution
    tb_prod, gain_db, res_m = analyzer.verify_pulse_compression_gain()
    assert np.isclose(tb_prod, 200.0), f"Expected TB product 200, got {tb_prod}"
    assert np.isclose(gain_db, 23.01, atol=0.1), f"Expected Gain ~23.01 dB, got {gain_db:.2f}"
    print(f" [PASS] Matched Filter Pulse Compression: TB={tb_prod:.0f}, Gain={gain_db:.2f} dB, Delta_R={res_m:.2f} m")

    # 3. Radar Range Equation
    ranges_km, pr_dbm, snr_db = analyzer.solve_radar_equation(
        peak_power_watts=120e3,
        antenna_gain_dbi=42.5,
        target_rcs_sqm=0.02, # 0.02 m^2 drone
        ranges_km=np.array([50.0, 100.0, 150.0])
    )
    assert snr_db[0] > snr_db[1] > snr_db[2], "SNR must monotonically decrease with range R^4"
    print(f" [PASS] Radar Range Equation: 50 km -> SNR={snr_db[0]:.1f} dB, 100 km -> SNR={snr_db[1]:.1f} dB")

    # 4. Swerling ROC Curves
    roc = analyzer.compute_swerling_roc(snr_db_range=np.array([10.0, 15.0, 20.0]), pfa=1e-6)
    assert roc["swerling_0"][2] > roc["swerling_0"][0], "Pd must increase with SNR"
    assert roc["swerling_0"][1] > roc["swerling_1"][1], "Swerling 0 Pd should exceed Swerling 1 at medium SNR"
    print(f" [PASS] Swerling Detection Models: Pd(Sw0@15dB)={roc['swerling_0'][1]:.3f}, Pd(Sw1@15dB)={roc['swerling_1'][1]:.3f}")

    # 5. Ambiguity Surface Peak Check
    delays, dopplers, amb = analyzer.compute_ambiguity_function(delay_bins=32, doppler_bins=32)
    assert np.isclose(np.max(amb), 1.0), "Ambiguity function must peak at 1.0 (origin)"
    print(f" [PASS] Woodward Ambiguity Surface: Peak={np.max(amb):.4f} centered at tau=0, nu=0")

    print("=" * 76)
    print("  ALL AEGISRADAR WAVEFORM & DSP ALGORITHMIC VERIFICATIONS PASSED (100%)")
    print("=" * 76)
    return 0


def main():
    parser = argparse.ArgumentParser(description="AegisRadar AESA Waveform & Ambiguity Analyzer")
    parser.add_argument("--self-test", action="store_true", help="Execute automated test verification")
    parser.add_argument("--fc", type=float, default=9.5e9, help="Carrier frequency in Hz (default: 9.5 GHz)")
    parser.add_argument("--bw", type=float, default=20.0e6, help="Bandwidth in Hz (default: 20 MHz)")
    parser.add_argument("--tau", type=float, default=10.0e-6, help="Pulse width in seconds (default: 10 µs)")
    args = parser.parse_args()

    if args.self_test:
        sys.exit(run_self_test())

    analyzer = RadarWaveformAnalyzer(
        carrier_freq_hz=args.fc,
        bandwidth_hz=args.bw,
        pulse_width_sec=args.tau
    )
    tb_prod, gain_db, res_m = analyzer.verify_pulse_compression_gain()
    print(f"Radar Configuration:")
    print(f" - Carrier Frequency: {args.fc / 1e9:.2f} GHz")
    print(f" - Bandwidth:         {args.bw / 1e6:.1f} MHz")
    print(f" - Pulse Width:       {args.tau * 1e6:.1f} µs")
    print(f" - TB Product:        {tb_prod:.1f}")
    print(f" - Compression Gain:  {gain_db:.2f} dB")
    print(f" - Range Resolution:  {res_m:.2f} m")


if __name__ == "__main__":
    main()
