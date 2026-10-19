#pragma once

#include <complex>
#include <vector>
#include <cmath>
#include <string>
#include <algorithm>
#include <numeric>
#include <stdexcept>

namespace aegis::dsp {

using ComplexF = std::complex<float>;
using ComplexD = std::complex<double>;

constexpr double SPEED_OF_LIGHT = 299792458.0; // m/s
constexpr double PI = 3.14159265358979323846;

struct RadarConfig {
    double carrierFrequencyHz{9.5e9};     // 9.5 GHz X-Band
    double bandwidthHz{20.0e6};            // 20 MHz Chirp bandwidth
    double pulseWidthSec{10.0e-6};         // 10 microseconds
    double samplingRateHz{50.0e6};         // 50 MHz ADC rate
    double prfHz{12500.0};                 // 12.5 kHz PRF
    double peakPowerWatts{120.0e3};        // 120 kW
    double antennaGainDbi{42.5};           // 42.5 dBi
    std::size_t numPulsesInCpi{32};        // 32 pulses per CPI
};

struct MonopulseAngles {
    double azimuthErrorDeg{0.0};
    double elevationErrorDeg{0.0};
    double sumChannelMagnitude{0.0};
    double azDiffMagnitude{0.0};
    double elDiffMagnitude{0.0};
    bool validLock{false};
};

struct DetectionBlip {
    std::size_t rangeBinIndex{0};
    std::size_t dopplerBinIndex{0};
    double estimatedRangeMeters{0.0};
    double radialVelocityMps{0.0};
    double snrDb{0.0};
    double powerLinear{0.0};
};

/**
 * High-performance Radix-2 Cooley-Tukey FFT implementation in C++17
 */
void fftRadix2(std::vector<ComplexD>& data, bool inverse = false);

/**
 * Fast FIR / Matched Filter Pulse Compression Engine
 */
class PulseCompressionEngine {
public:
    explicit PulseCompressionEngine(const RadarConfig& cfg);

    [[nodiscard]] const std::vector<ComplexD>& getReferenceChirp() const noexcept {
        return referenceChirp_;
    }

    [[nodiscard]] double getTimeBandwidthProduct() const noexcept {
        return cfg_.bandwidthHz * cfg_.pulseWidthSec;
    }

    [[nodiscard]] double getRangeResolutionMeters() const noexcept {
        return SPEED_OF_LIGHT / (2.0 * cfg_.bandwidthHz);
    }

    /**
     * Executes matched filter pulse compression: Y(f) = X(f) * H*(f)
     */
    [[nodiscard]] std::vector<ComplexD> compress(const std::vector<ComplexD>& rxSignal) const;

private:
    RadarConfig cfg_;
    std::vector<ComplexD> referenceChirp_;
    void synthesizeChirp();
};

/**
 * 2D Pulse-Doppler & Moving Target Detector (MTD) Processor
 */
class PulseDopplerProcessor {
public:
    explicit PulseDopplerProcessor(const RadarConfig& cfg);

    /**
     * Applies 3-Pulse Clutter Delay-Line Canceller (MTI)
     * Transfer function: H(z) = (1 - z^-1)^2 = 1 - 2z^-1 + z^-2
     */
    [[nodiscard]] std::vector<ComplexD> applyThreePulseMti(const std::vector<ComplexD>& slowTimePulses) const;

    /**
     * Transforms 2D slow-time vs fast-time CPI matrix into Range-Doppler Matrix
     */
    [[nodiscard]] std::vector<std::vector<double>> computeRangeDopplerMap(
        const std::vector<std::vector<ComplexD>>& cpiMatrix) const;

    [[nodiscard]] double getBlindSpeedMps() const noexcept {
        const double lambda = SPEED_OF_LIGHT / cfg_.carrierFrequencyHz;
        return (lambda * cfg_.prfHz) / 2.0;
    }

private:
    RadarConfig cfg_;
};

/**
 * Amplitude-Comparison Monopulse Angle Tracking Sub-aperture Processor
 */
class MonopulseTracker {
public:
    MonopulseTracker(double azSlope = 1.6, double elSlope = 1.6) noexcept
        : azSlope_(azSlope), elSlope_(elSlope) {}

    /**
     * Computes Sum (Sigma), Azimuth Diff (Delta_az), and Elevation Diff (Delta_el)
     * from 4 AESA array quadrant sub-apertures:
     * Quadrant layout:
     *   [A] [B]
     *   [C] [D]
     */
    [[nodiscard]] MonopulseAngles processQuadrants(
        ComplexD qA, ComplexD qB, ComplexD qC, ComplexD qD) const;

private:
    double azSlope_;
    double elSlope_;
};

/**
 * 2D Cell-Averaging Constant False Alarm Rate (CA-CFAR) Detector
 */
class CaCfarDetector {
public:
    CaCfarDetector(std::size_t numGuard = 2, std::size_t numRef = 4, double pfa = 1e-6)
        : numGuard_(numGuard), numRef_(numRef), pfa_(pfa) {}

    [[nodiscard]] std::vector<DetectionBlip> detect2D(
        const std::vector<std::vector<double>>& rangeDopplerMap,
        const RadarConfig& cfg) const;

private:
    std::size_t numGuard_;
    std::size_t numRef_;
    double pfa_;
};

} // namespace aegis::dsp
