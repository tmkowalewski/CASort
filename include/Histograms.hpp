#ifndef HISTOGRAMS_HPP
#define HISTOGRAMS_HPP

// C++ Includes
#include <memory>

// ROOT Includes
#include <TH1D.h>
#include <TH2D.h>
#include <ROOT/TThreadedObject.hxx>

// Project Includes
#include "Configuration.hpp"

template <typename T>
struct Histogram
{
    template <typename... Args>
    explicit Histogram(Args &&...args)
        : fHistogram(std::make_unique<ROOT::TThreadedObject<T>>(std::forward<Args>(args)...)) {}

    auto GetThreadLocalPtr() { return fHistogram->Get(); }

    auto Merge() { return fHistogram->Merge(); }

    auto Write() { return fHistogram->Merge()->Write(); }

    std::unique_ptr<ROOT::TThreadedObject<T>> fHistogram;
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
#if PROCESS_CLOVER_CROSS

    // Raw Histos
    // Raw Histos
    auto cc_amp = Histogram<TH2D>("cc_amp", "Clover Cross Amplitude (Raw Data);ADC;Channel;Counts/Bin", kDigitizerBins, 0, kDigitizerBins, kDigitizerChannels, 0, kDigitizerChannels);
    auto cc_cht = Histogram<TH2D>("cc_cht", "Clover Cross Channel Time (Raw Data);ADC;Channel;Counts/Bin", kDigitizerBins, 0, (kDigitizerBins)*kNsPerBin, kDigitizerChannels, 0, kDigitizerChannels);
    auto cc_plu = Histogram<TH2D>("cc_plu", "Clover Cross Pile-Up;Pile-Up Multiplicity;Channel;Counts/Bin", kDigitizerBins, 0, kDigitizerBins, kDigitizerChannels, 0, kDigitizerChannels);
    auto cc_mdt = Histogram<TH1D>("cc_mdt", "Clover Cross Module Time;Time (ns);Counts/Bin", kDigitizerBins, 0, (kDigitizerBins)*kNsPerBin);
    auto cc_trt = Histogram<TH2D>("cc_trt", "Clover Cross Trigger Time;Time (ns);Trigger ID;Counts/Bin", kDigitizerBins, 0, (kDigitizerBins)*kNsPerBin, 2, 0, 2);

    // Calibrated Histos
    auto cc_E = Histogram<TH2D>("cc_E", "Clover Cross Energy;Energy (keV);Channel;Counts/Bin", kMaxEnergy / kEnergyPerBin, 0, kMaxEnergy, kDigitizerChannels, 0, kDigitizerChannels);

    // Sum Histos
    auto cc_sum = Histogram<TH2D>("cc_sum", "Clover Cross Energy (Detector Sum);Energy (keV);Detector;Counts/Bin", kMaxEnergy / kEnergyPerBin, 0, kMaxEnergy, kDigitizerChannels / 4, 0, kDigitizerChannels / 4);

    // Addback Histos
    auto cc_abE = Histogram<TH2D>("cc_abE", "Clover Cross Energy (Detector Addback);Energy (keV);Detector;Counts/Bin", kMaxEnergy / kEnergyPerBin, 0, kMaxEnergy, kDigitizerChannels / 4, 0, kDigitizerChannels / 4);
    auto cc_abM = Histogram<TH1D>("cc_abM", "Clover Cross Addback Multiplicity;Multiplicity;Counts/Bin", 4, 0, 4);

#endif // PROCESS_CLOVER_CROSS
/* #endregion */

/* #region clover_back */
#if PROCESS_CLOVER_BACK

    // Raw Histos
    auto clover_back_amp_raw = Histogram<TH2D>("clover_back_amp_raw", "Clover Back Amplitude (Raw Data);ADC;Channel;Counts/Bin", kDigitizerBins, 0, kDigitizerBins, kDigitizerChannels, 0, kDigitizerChannels);
    auto clover_back_cht_raw = Histogram<TH2D>("clover_back_cht_raw", "Clover Back Channel Time (Raw Data);ADC;Channel;Counts/Bin", kDigitizerBins, 0, kDigitizerBins, kDigitizerChannels, 0, kDigitizerChannels);
    auto clover_back_plu = Histogram<TH2D>("clover_back_plu", "Clover Back Pileup (Raw Data);Multiplicity;Channel;Counts/Bin", kDigitizerBins, 0, kDigitizerBins, kDigitizerChannels, 0, kDigitizerChannels);

    auto clover_back_mdt_raw = Histogram<TH1D>("clover_back_mdt_raw", "Clover Back Module Time (Raw Data);ADC;Counts/Bin", kDigitizerBins, 0, kDigitizerBins);
    auto clover_back_trt_raw = Histogram<TH2D>("clover_back_trt_raw", "Clover Back Trigger Time (Raw Data);ADC;Trigger ID;Counts/Bin", kDigitizerBins, 0, kDigitizerBins, 2, 0, 2);

    // Calibrated Histos
    auto clover_back_E = Histogram<TH2D>("clover_back_E", "Clover Back Energy;Energy (keV);Channel;Counts/Bin", kMaxEnergy / kEnergyPerBin, 0, kMaxEnergy, kDigitizerChannels, 0, kDigitizerChannels);
    auto clover_back_cht = Histogram<TH2D>("clover_back_cht", "Clover Back Channel Time;Channel Time (ns);Counts/Bin", kDigitizerBins, 0, (kDigitizerBins)*kNsPerBin, kDigitizerChannels, 0, kDigitizerChannels);
    auto clover_back_mdt = Histogram<TH1D>("clover_back_mdt", "Clover Back Module Time; Module Time (ns); Counts/Bin", kDigitizerBins, 0, (kDigitizerBins)*kNsPerBin);
    auto clover_back_trt = Histogram<TH2D>("clover_back_trt", "Clover Back Trigger Time; Trigger Time (ns); Trigger ID; Counts/Bin", kDigitizerBins, 0, (kDigitizerBins)*kNsPerBin, 2, 0, 2);

    // Sum Histos
    auto clover_back_sum = Histogram<TH2D>("clover_back_sum", "Clover Back Energy (Detector Sum);Energy (keV);Detector;Counts/Bin", kMaxEnergy / kEnergyPerBin, 0, kMaxEnergy, kDigitizerChannels / 4, 0, kDigitizerChannels / 4);

    // Addback Histos
    auto clover_back_addback = Histogram<TH2D>("clover_back_addback", "Clover Back Energy (Detector Addback);Energy (keV);Detector;Counts/Bin", kMaxEnergy / kEnergyPerBin, 0, kMaxEnergy, kDigitizerChannels / 4, 0, kDigitizerChannels / 4);
    auto clover_back_addback_mult = Histogram<TH1D>("clover_back_addback_mult", "Clover Back Addback Multiplicity;Multiplicity;Counts/Bin", 4, 0, 4);

#endif // PROCESS_CLOVER_BACK
/* #endregion */

/* #region pos_sig */
#if PROCESS_POS_SIG

    // Raw Histos
    auto pos_sig_amp_raw = Histogram<TH2D>("pos_sig_amp_raw", "Positive Signal Amplitude (Raw Data);ADC;Channel;Counts/Bin", kDigitizerBins, 0, kDigitizerBins, kDigitizerChannels, 0, kDigitizerChannels);
    auto pos_sig_cht_raw = Histogram<TH2D>("pos_sig_cht_raw", "Positive Signal Channel Time (Raw Data);ADC;Channel;Counts/Bin", kDigitizerBins, 0, kDigitizerBins, kDigitizerChannels, 0, kDigitizerChannels);
    auto pos_sig_plu = Histogram<TH2D>("pos_sig_plu", "Positive Signal Pileup (Raw Data);ADC;Channel;Counts/Bin", kDigitizerBins, 0, kDigitizerBins, kDigitizerChannels, 0, kDigitizerChannels);

    auto pos_sig_mdt_raw = Histogram<TH1D>("pos_sig_mdt_raw", "Positive Signal Module Time (Raw Data);ADC;Counts/Bin", kDigitizerBins, 0, kDigitizerBins);
    auto pos_sig_trt_raw = Histogram<TH2D>("pos_sig_trt_raw", "Positive Signal Trigger Time (Raw Data);ADC;Trigger ID;Counts/Bin", kDigitizerBins, 0, kDigitizerBins, 2, 0, 2);

    // Calibrated Histos
    auto pos_sig_E = Histogram<TH2D>("pos_sig_E", "Pos Sig Energy;Energy (keV),Crystal,Counts/Bin", kMaxEnergy / kEnergyPerBin, 0, kMaxEnergy, kDigitizerChannels, 0, kDigitizerChannels);
    auto pos_sig_cht = Histogram<TH2D>("pos_sig_cht", "Pos Sig Channel Time;Channel Time (ns);Channel;Counts/Bin", kDigitizerBins, 0, (kDigitizerBins)*kNsPerBin, kDigitizerChannels, 0, kDigitizerChannels);
    auto pos_sig_mdt = Histogram<TH1D>("pos_sig_mdt", "Positive Signal Module Time;Module Time (ns);Counts/Bin", kDigitizerBins, 0, (kDigitizerBins)*kNsPerBin);
    auto pos_sig_trt = Histogram<TH2D>("pos_sig_trt", "Positive Signal Trigger Time;Trigger Time (ns);Trigger ID;Counts/Bin", kDigitizerBins, 0, (kDigitizerBins)*kNsPerBin, 2, 0, 2);

    // Sum Histos
    auto pos_sig_sum = Histogram<TH2D>("pos_sig_sum", "Pos Sig Energy (Detector Sum);Energy (keV);Detector;Counts/bin", kMaxEnergy / kEnergyPerBin, 0, kMaxEnergy, kDigitizerChannels / 4, 0, kDigitizerChannels / 4);

    // Addback Histos
    auto pos_sig_addback = Histogram<TH2D>("pos_sig_addback", "Pos Sig Energy (Detctor Addback);Energy (keV);Detetcor;Counts/Bin", kMaxEnergy / kEnergyPerBin, 0, kMaxEnergy, kDigitizerChannels / 4, 0, kDigitizerChannels / 4);
    auto pos_sig_addback_mult = Histogram<TH1D>("pos_sig_addback_mult", "Pos Sig Addback Multiplicity;Multiplicity;Counts/Bin", 4, 0, 4);

#endif // PROCESS_POS_SIG

/* #region cebr_all */
#if PROCESS_CEBR_ALL

    // Raw Histos
    auto cebr_all_inl_raw = Histogram<TH2D>("cebr_all_inl_raw", "CeBr All Integration Long (Raw Data);ADC;Channel;Counts/Bin", kDigitizerBins, 0, kDigitizerBins, kDigitizerChannels, 0, kDigitizerChannels);
    auto cebr_all_cht_raw = Histogram<TH2D>("cebr_all_cht_raw", "CeBr All Channel Time (Raw Data);ADC;Channel;Counts/Bin", kDigitizerBins, 0, kDigitizerBins, kDigitizerChannels, 0, kDigitizerChannels);
    auto cebr_all_ins_raw = Histogram<TH2D>("cebr_all_ins_raw", "CeBr All Integration Short (Raw Data);ADC;Channel;Counts/Bin", kDigitizerBins, 0, kDigitizerBins, kDigitizerChannels, 0, kDigitizerChannels);

    auto cebr_all_mdt_raw = Histogram<TH1D>("cebr_all_mdt_raw", "CeBr All Module Time (Raw Data);ADC;Counts/Bin", kDigitizerBins, 0, kDigitizerBins);
    auto cebr_all_trt_raw = Histogram<TH2D>("cebr_all_trt_raw", "CeBr All Trigger Time (Raw Data);ADC;Trigger ID;Counts/Bin", kDigitizerBins, 0, kDigitizerBins, 2, 0, 2);

    // Calibrated Histos
    auto cebr_all_El = Histogram<TH2D>("cebr_all_El", "CeBr All Energy (Long Integration);Energy (keV);Channel;Counts/Bin", kMaxEnergy / kEnergyPerBin, 0, kMaxEnergy, kDigitizerChannels, 0, kDigitizerChannels);
    auto cebr_all_cht = Histogram<TH2D>("cebr_all_cht", "CeBr All Channel Time; Channel Time (ns);Channel;Counts/Bin", kDigitizerBins, 0, (kDigitizerBins)*kNsPerBin, kDigitizerChannels, 0, kDigitizerChannels);
    auto cebr_all_Es = Histogram<TH2D>("cebr_all_Es", "CeBr All Energy (Short Integration); Energy (keV);Channel;Counts/Bin", kDigitizerBins, 0, kDigitizerBins, kDigitizerChannels, 0, kDigitizerChannels);

    auto cebr_all_mdt = Histogram<TH1D>("cebr_all_mdt", "CeBr All Module Time;Module Time (ns);Counts/Bin", kDigitizerBins, 0, (kDigitizerBins)*kNsPerBin);
    auto cebr_all_trt = Histogram<TH2D>("cebr_all_trt", "CeBr All Trigger Time;Trigger Time (ns);Trigger ID;Counts/Bin", kDigitizerBins, 0, (kDigitizerBins)*kNsPerBin, 2, 0, 2);

#endif // PROCESS_CEBR_ALL
    /* #endregion */

} // namespace Histograms

#endif // HISTOGRAMS_HPP