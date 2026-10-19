#include "aesa_pulse_doppler.hpp"
#include <iostream>

namespace aegis::dsp {

static std::size_t nextPowerOfTwo(std::size_t n) {
    if (n <= 1) return 1;
    std::size_t p = 1;
    while (p < n) {
        p <<= 1;
    }
    return p;
}

void fftRadix2(std::vector<ComplexD>& data, bool inverse) {
    const std::size_t n = data.size();
    if ((n & (n - 1)) != 0) {
        throw std::invalid_argument("FFT length must be an exact power of two");
    }

    // Bit-reversal permutation
    std::size_t j = 0;
    for (std::size_t i = 0; i < n - 1; ++i) {
        if (i < j) {
            std::swap(data[i], data[j]);
        }
        std::size_t k = n >> 1;
        while (k <= j) {
            j -= k;
            k >>= 1;
        }
        j += k;
    }

    // Butterfly computations
    const double sign = inverse ? 1.0 : -1.0;
    for (std::size_t len = 2; len <= n; len <<= 1) {
        const std::size_t halfLen = len >> 1;
        const double angle = (sign * 2.0 * PI) / static_cast<double>(len);
        const ComplexD wBase(std::cos(angle), std::sin(angle));

        for (std::size_t i = 0; i < n; i += len) {
            ComplexD w(1.0, 0.0);
            for (std::size_t k = 0; k < halfLen; ++k) {
                const ComplexD u = data[i + k];
                const ComplexD v = data[i + k + halfLen] * w;
                data[i + k] = u + v;
                data[i + k + halfLen] = u - v;
                w *= wBase;
            }
        }
    }

    // Scale for inverse transform
    if (inverse) {
        const double invN = 1.0 / static_cast<double>(n);
        for (std::size_t i = 0; i < n; ++i) {
            data[i] *= invN;
        }
    }
}

PulseCompressionEngine::PulseCompressionEngine(const RadarConfig& cfg)
    : cfg_(cfg) {
    synthesizeChirp();
}

void PulseCompressionEngine::synthesizeChirp() {
    const std::size_t numSamples = static_cast<std::size_t>(
        std::floor(cfg_.pulseWidthSec * cfg_.samplingRateHz));
    referenceChirp_.resize(numSamples);

    const double dt = 1.0 / cfg_.samplingRateHz;
    const double chirpRate = cfg_.bandwidthHz / cfg_.pulseWidthSec;

    for (std::size_t i = 0; i < numSamples; ++i) {
        const double t = -cfg_.pulseWidthSec / 2.0 + static_cast<double>(i) * dt;
        const double phase = PI * chirpRate * t * t;
        referenceChirp_[i] = ComplexD(std::cos(phase), std::sin(phase));
    }
}

std::vector<ComplexD> PulseCompressionEngine::compress(
    const std::vector<ComplexD>& rxSignal) const {
    const std::size_t minLen = rxSignal.size() + referenceChirp_.size() - 1;
    const std::size_t fftLen = nextPowerOfTwo(minLen);

    std::vector<ComplexD> paddedRx(fftLen, ComplexD(0.0, 0.0));
    std::copy(rxSignal.begin(), rxSignal.end(), paddedRx.begin());

    std::vector<ComplexD> paddedRef(fftLen, ComplexD(0.0, 0.0));
    std::copy(referenceChirp_.begin(), referenceChirp_.end(), paddedRef.begin());

    fftRadix2(paddedRx, false);
    fftRadix2(paddedRef, false);

    // Fast frequency-domain matched filter multiplication
    for (std::size_t i = 0; i < fftLen; ++i) {
        paddedRx[i] = paddedRx[i] * std::conj(paddedRef[i]);
    }

    fftRadix2(paddedRx, true);
    paddedRx.resize(rxSignal.size());
    return paddedRx;
}

PulseDopplerProcessor::PulseDopplerProcessor(const RadarConfig& cfg)
    : cfg_(cfg) {}

std::vector<ComplexD> PulseDopplerProcessor::applyThreePulseMti(
    const std::vector<ComplexD>& slowTimePulses) const {
    if (slowTimePulses.size() < 3) {
        return slowTimePulses;
    }

    std::vector<ComplexD> output;
    output.reserve(slowTimePulses.size() - 2);

    for (std::size_t i = 2; i < slowTimePulses.size(); ++i) {
        // y[n] = x[n] - 2*x[n-1] + x[n-2]
        const ComplexD y = slowTimePulses[i] - 2.0 * slowTimePulses[i - 1] + slowTimePulses[i - 2];
        output.push_back(y);
    }
    return output;
}

std::vector<std::vector<double>> PulseDopplerProcessor::computeRangeDopplerMap(
    const std::vector<std::vector<ComplexD>>& cpiMatrix) const {
    if (cpiMatrix.empty()) return {};

    const std::size_t numPulses = cpiMatrix.size();
    const std::size_t numRangeBins = cpiMatrix[0].size();
    const std::size_t dopplerFftSize = nextPowerOfTwo(std::max(numPulses, std::size_t{16}));

    // Output matrix: [numRangeBins x dopplerFftSize]
    std::vector<std::vector<double>> rdMap(
        numRangeBins, std::vector<double>(dopplerFftSize, 0.0));

    // Process each range gate across slow-time pulses
    std::vector<ComplexD> slowTime(dopplerFftSize, ComplexD(0.0, 0.0));

    for (std::size_t r = 0; r < numRangeBins; ++r) {
        // Apply Hamming window across pulses
        for (std::size_t p = 0; p < dopplerFftSize; ++p) {
            if (p < numPulses) {
                const double w = 0.54 - 0.46 * std::cos((2.0 * PI * p) / static_cast<double>(numPulses - 1));
                slowTime[p] = cpiMatrix[p][r] * w;
            } else {
                slowTime[p] = ComplexD(0.0, 0.0);
            }
        }

        fftRadix2(slowTime, false);

        // FFT shift & power conversion
        const std::size_t half = dopplerFftSize / 2;
        for (std::size_t d = 0; d < dopplerFftSize; ++d) {
            const std::size_t shiftedIdx = (d + half) % dopplerFftSize;
            const double magSq = std::norm(slowTime[shiftedIdx]);
            rdMap[r][d] = 10.0 * std::log10(std::max(1e-12, magSq));
        }
    }

    return rdMap;
}

MonopulseAngles MonopulseTracker::processQuadrants(
    ComplexD qA, ComplexD qB, ComplexD qC, ComplexD qD) const {
    // Sum channel Sigma: A + B + C + D
    const ComplexD sigma = qA + qB + qC + qD;
    // Azimuth diff Delta_az: (A + C) - (B + D)
    const ComplexD deltaAz = (qA + qC) - (qB + qD);
    // Elevation diff Delta_el: (A + B) - (C + D)
    const ComplexD deltaEl = (qA + qB) - (qC + qD);

    const double sigmaMag = std::abs(sigma);
    MonopulseAngles res;
    res.sumChannelMagnitude = sigmaMag;
    res.azDiffMagnitude = std::abs(deltaAz);
    res.elDiffMagnitude = std::abs(deltaEl);

    if (sigmaMag < 1e-9) {
        res.validLock = false;
        return res;
    }

    // Monopulse ratio real parts scaled by beam monopulse slope
    const ComplexD ratioAz = deltaAz / sigma;
    const ComplexD ratioEl = deltaEl / sigma;

    res.azimuthErrorDeg = ratioAz.real() * azSlope_;
    res.elevationErrorDeg = ratioEl.real() * elSlope_;
    res.validLock = true;
    return res;
}

std::vector<DetectionBlip> CaCfarDetector::detect2D(
    const std::vector<std::vector<double>>& rangeDopplerMap,
    const RadarConfig& cfg) const {
    std::vector<DetectionBlip> detections;
    if (rangeDopplerMap.empty()) return detections;

    const std::size_t numR = rangeDopplerMap.size();
    const std::size_t numD = rangeDopplerMap[0].size();

    const double dtFast = 1.0 / cfg.samplingRateHz;
    const double rangeRes = (SPEED_OF_LIGHT * dtFast) / 2.0;
    const double lambda = SPEED_OF_LIGHT / cfg.carrierFrequencyHz;
    const double prf = cfg.prfHz;

    // CA-CFAR alpha factor
    const std::size_t totalRefCells = 2 * numRef_;
    const double alpha = static_cast<double>(totalRefCells) *
        (std::pow(pfa_, -1.0 / static_cast<double>(totalRefCells)) - 1.0);

    for (std::size_t r = numGuard_ + numRef_; r < numR - (numGuard_ + numRef_); ++r) {
        for (std::size_t d = 0; d < numD; ++d) {
            double refSumLinear = 0.0;
            std::size_t refCount = 0;

            // Average training cells along range dimension
            for (std::size_t ir = r - (numGuard_ + numRef_); ir <= r + numGuard_ + numRef_; ++ir) {
                if (ir >= r - numGuard_ && ir <= r + numGuard_) {
                    continue; // Skip guard cells and CUT
                }
                const double cellPwrLinear = std::pow(10.0, rangeDopplerMap[ir][d] / 10.0);
                refSumLinear += cellPwrLinear;
                refCount++;
            }

            if (refCount == 0) continue;

            const double noiseEstimate = refSumLinear / static_cast<double>(refCount);
            const double adaptiveThreshold = noiseEstimate * alpha;
            const double cutPwr = std::pow(10.0, rangeDopplerMap[r][d] / 10.0);

            if (cutPwr > adaptiveThreshold) {
                DetectionBlip blip;
                blip.rangeBinIndex = r;
                blip.dopplerBinIndex = d;
                blip.estimatedRangeMeters = static_cast<double>(r) * rangeRes;

                // Doppler to radial velocity: fd = (d - numD/2) * (PRF / numD)
                const double fd = (static_cast<double>(d) - static_cast<double>(numD) / 2.0) *
                                  (prf / static_cast<double>(numD));
                blip.radialVelocityMps = (fd * lambda) / 2.0;
                blip.powerLinear = cutPwr;
                blip.snrDb = 10.0 * std::log10(std::max(1e-3, cutPwr / noiseEstimate));

                detections.push_back(blip);
            }
        }
    }

    return detections;
}

} // namespace aegis::dsp
