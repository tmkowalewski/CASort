#ifndef HISTOGRAMS_HPP
#define HISTOGRAMS_HPP

// C++ Includes
#include <memory>

// ROOT Includes
#include <TH1D.h>
#include <TH2D.h>
#include <ROOT/TThreadedObject.hxx>

template <typename T>
struct Histogram
{
    template <typename... Args>
    explicit Histogram(Args &&...args)
        : fHistogram(std::make_shared<ROOT::TThreadedObject<T>>(std::forward<Args>(args)...)) {}

    auto MakePtr() { return fHistogram->Get(); }

    auto Merge() { return fHistogram->Merge(); }

    auto Write() { return fHistogram->Merge()->Write(); }

    std::shared_ptr<ROOT::TThreadedObject<T>> fHistogram;
};

namespace Histograms
{
    // Constants
    static const constexpr double kMaxEnergy = 10000.0;  // Maximum energy for histograms in keV
    static const constexpr double kEnergyPerBin = 0.25;  // Energy per bin in keV
    static const constexpr double kNsPerBin = 0.098;     // Conversion factor from bin to nanoseconds
    static const constexpr int kDigitizerBins = 1 << 16; // Number of bins in the digitizer (16-bit)
    static const constexpr int kDigitizerChannels = 16;  // Number of channels in digitizer

    /* #region clover_cross */

    // Raw Histos
    auto clover_cross_amp_raw = Histogram<TH2D>("clover_cross_amp_raw", "Clover Cross Amplitude (Raw Data);ADC;Channel;Counts/Bin", kDigitizerBins, 0, kDigitizerBins, kDigitizerChannels, 0, kDigitizerChannels);

    /* #endregion */

}

#endif // HISTOGRAMS_HPP