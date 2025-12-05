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
    // Raw Histos
    auto clover_cross_amp_raw = Histogram<TH2D>("clover_cross_amp_raw", "Clover Cross Amplitude (Raw Data);ADC;Channel;Counts/Bin", kDigitizerBins, 0, kDigitizerBins, kDigitizerChannels, 0, kDigitizerChannels);
    auto clover_cross_cht_raw = Histogram<TH2D>("clover_cross_cht_raw", "Clover Cross Channel Time (Raw Data);ADC;Channel;Counts/Bin", kDigitizerBins, 0, kDigitizerBins, kDigitizerChannels, 0, kDigitizerChannels);
    auto clover_cross_plu = Histogram<TH2D>("clover_cross_plu", "Clover Cross Pileup (Raw Data);Multiplicity;Channel;Counts/Bin", kDigitizerBins, 0, kDigitizerBins, kDigitizerChannels, 0, kDigitizerChannels);

    auto clover_cross_mdt_raw = Histogram<TH1D>("clover_cross_mdt_raw", "Clover Cross Module Time (Raw Data);ADC;Counts/Bin", kDigitizerBins, 0, kDigitizerBins);
    auto clover_cross_trt_raw = Histogram<TH2D>("clover_cross_trt_raw", "Clover Cross Trigger Time (Raw Data);ADC;Trigger ID;Counts/Bin", kDigitizerBins, 0, kDigitizerBins, 2, 0, 2);

    // Calibrated Histos
    auto clover_cross_E = Histogram<TH2D>("clover_cross_E", "Clover Cross Energy;Energy (keV);Channel;Counts/Bin", kMaxEnergy / kEnergyPerBin, 0, kMaxEnergy, kDigitizerChannels, 0, kDigitizerChannels);
    auto clover_cross_cht = Histogram<TH2D>("clover_cross_cht", "Clover Cross Channel Time;Channel Time (ns);Counts/Bin", kDigitizerBins, 0, (kDigitizerBins)*kNsPerBin, kDigitizerChannels, 0, kDigitizerChannels);
    auto clover_cross_mdt = Histogram<TH1D>("clover_cross_mdt", "Clover Cross Module Time; Module Time (ns); Counts/Bin", kDigitizerBins, 0, (kDigitizerBins)*kNsPerBin);
    auto clover_cross_trt = Histogram<TH2D>("clover_cross_trt", "Clover Cross Trigger Time; Trigger Time (ns); Trigger ID; Counts/Bin", kDigitizerBins, 0, (kDigitizerBins)*kNsPerBin, 2, 0, 2);

    // Sum Histos
    auto clover_cross_sum = Histogram<TH2D>("clover_cross_sum", "Clover Cross Energy (Detector Sum);Energy (keV);Detector;Counts/Bin", kMaxEnergy / kEnergyPerBin, 0, kMaxEnergy, kDigitizerChannels / 4, 0, kDigitizerChannels / 4);

    // Addback Histos
    auto clover_cross_addback = Histogram<TH2D>("clover_cross_addback", "Clover Cross Energy (Detector Addback);Energy (keV);Detector;Counts/Bin", kMaxEnergy / kEnergyPerBin, 0, kMaxEnergy, kDigitizerChannels / 4, 0, kDigitizerChannels / 4);
    auto clover_cross_addback_mult = Histogram<TH1D>("clover_cross_addback_mult", "Clover Cross Addback Multiplicity;Multiplicity;Counts/Bin", 4, 0, 4);

    /* #endregion */

}

#endif // HISTOGRAMS_HPP