#ifndef CAHISTOGRAMS_HPP
#define CAHISTOGRAMS_HPP

// C++ Includes
#include <memory>

// ROOT Includes
#include <TH1D.h>
#include <TH2D.h>
#include <ROOT/TThreadedObject.hxx>

// Project Includes
#include "CAConfiguration.hpp"
#include "TCAHistogram.hpp"

namespace CAHistograms
{
    // Constants
    inline constexpr double kMaxEnergy = 10000.0;           // Maximum energy for histograms in keV
    inline constexpr double kEnergyPerBin = 0.25;           // Energy per bin in keV
    inline constexpr double kNsPerBin = 0.098;              // Conversion factor from bin to nanoseconds
    inline constexpr unsigned int kDigitizerBins = 1 << 16; // Number of bins in the digitizer (16-bit)
    inline constexpr unsigned int kDigitizerChannels = 16;  // Number of channels in digitizer

    #if PROCESS_CLOVER_CROSS

    // Raw Hists
    auto cc_amp = TCAHistogram<TH2D>("cc_amp", "Clover Cross Amplitude (Raw Data);ADC;Channel;Counts/Bin", kDigitizerBins, 0, kDigitizerBins, kDigitizerChannels, 0, kDigitizerChannels);
    auto cc_cht = TCAHistogram<TH2D>("cc_cht", "Clover Cross Channel Time (Raw Data);ADC;Channel;Counts/Bin", kDigitizerBins, 0, (kDigitizerBins)*kNsPerBin, kDigitizerChannels, 0, kDigitizerChannels);
    auto cc_plu = TCAHistogram<TH2D>("cc_plu", "Clover Cross Pile-Up;Pile-Up Multiplicity;Channel;Counts/Bin", kDigitizerBins, 0, kDigitizerBins, kDigitizerChannels, 0, kDigitizerChannels);
    auto cc_mdt = TCAHistogram<TH1D>("cc_mdt", "Clover Cross Module Time;Time (ns);Counts/Bin", kDigitizerBins, 0, (kDigitizerBins)*kNsPerBin);
    auto cc_trt = TCAHistogram<TH2D>("cc_trt", "Clover Cross Trigger Time;Time (ns);Trigger ID;Counts/Bin", kDigitizerBins, 0, (kDigitizerBins)*kNsPerBin, 2, 0, 2);

    // Calibrated Hists
    auto cc_xtE = TCAHistogram<TH2D>("cc_xtE", "Clover Cross Energy;Energy (keV);Channel;Counts/Bin", kMaxEnergy / kEnergyPerBin, 0, kMaxEnergy, kDigitizerChannels, 0, kDigitizerChannels);

    // Sum Hists
    auto cc_sum = TCAHistogram<TH2D>("cc_sum", "Clover Cross Energy (Detector Sum);Energy (keV);Detector;Counts/Bin", kMaxEnergy / kEnergyPerBin, 0, kMaxEnergy, kDigitizerChannels / 4, 0, kDigitizerChannels / 4);

    // Addback Hists
    auto cc_abE = TCAHistogram<TH2D>("cc_abE", "Clover Cross Energy (Detector Addback);Energy (keV);Detector;Counts/Bin", kMaxEnergy / kEnergyPerBin, 0, kMaxEnergy, kDigitizerChannels / 4, 0, kDigitizerChannels / 4);
    auto cc_abM = TCAHistogram<TH1D>("cc_abM", "Clover Cross Addback Multiplicity;Multiplicity;Counts/Bin", 4, 0, 4);

    #endif // PROCESS_CLOVER_CROSS

    #if PROCESS_CLOVER_BACK

    // Raw Hists
    auto cb_amp = TCAHistogram<TH2D>("cb_amp", "Clover Back Amplitude (Raw Data);ADC;Channel;Counts/Bin", kDigitizerBins, 0, kDigitizerBins, kDigitizerChannels, 0, kDigitizerChannels);
    auto cb_cht = TCAHistogram<TH2D>("cb_cht", "Clover Back Channel Time (Raw Data);ADC;Channel;Counts/Bin", kDigitizerBins, 0, (kDigitizerBins)*kNsPerBin, kDigitizerChannels, 0, kDigitizerChannels);
    auto cb_plu = TCAHistogram<TH2D>("cb_plu", "Clover Back Pile-Up;Pile-Up Multiplicity;Channel;Counts/Bin", kDigitizerBins, 0, kDigitizerBins, kDigitizerChannels, 0, kDigitizerChannels);
    auto cb_mdt = TCAHistogram<TH1D>("cb_mdt", "Clover Back Module Time;Time (ns);Counts/Bin", kDigitizerBins, 0, (kDigitizerBins)*kNsPerBin);
    auto cb_trt = TCAHistogram<TH2D>("cb_trt", "Clover Back Trigger Time;Time (ns);Trigger ID;Counts/Bin", kDigitizerBins, 0, (kDigitizerBins)*kNsPerBin, 2, 0, 2);

    // Calibrated Hists
    auto cb_xtE = TCAHistogram<TH2D>("cb_xtE", "Clover Back Energy;Energy (keV);Channel;Counts/Bin", kMaxEnergy / kEnergyPerBin, 0, kMaxEnergy, kDigitizerChannels, 0, kDigitizerChannels);

    // Sum Hists
    auto cb_sum = TCAHistogram<TH2D>("cb_sum", "Clover Back Energy (Detector Sum);Energy (keV);Detector;Counts/Bin", kMaxEnergy / kEnergyPerBin, 0, kMaxEnergy, kDigitizerChannels / 4, 0, kDigitizerChannels / 4);

    // Addback Hists
    auto cb_abE = TCAHistogram<TH2D>("cb_abE", "Clover Back Energy (Detector Addback);Energy (keV);Detector;Counts/Bin", kMaxEnergy / kEnergyPerBin, 0, kMaxEnergy, kDigitizerChannels / 4, 0, kDigitizerChannels / 4);
    auto cb_abM = TCAHistogram<TH1D>("cb_abM", "Clover Back Addback Multiplicity;Multiplicity;Counts/Bin", 4, 0, 4);


    #endif // PROCESS_CLOVER_BACK

    #if PROCESS_POS_SIG

    // Raw Hists
    auto pos_sig_amp_raw = TCAHistogram<TH2D>("pos_sig_amp_raw", "Positive Signal Amplitude (Raw Data);ADC;Channel;Counts/Bin", kDigitizerBins, 0, kDigitizerBins, kDigitizerChannels, 0, kDigitizerChannels);
    auto pos_sig_cht_raw = TCAHistogram<TH2D>("pos_sig_cht_raw", "Positive Signal Channel Time (Raw Data);ADC;Channel;Counts/Bin", kDigitizerBins, 0, kDigitizerBins, kDigitizerChannels, 0, kDigitizerChannels);
    auto pos_sig_plu = TCAHistogram<TH2D>("pos_sig_plu", "Positive Signal Pileup (Raw Data);ADC;Channel;Counts/Bin", kDigitizerBins, 0, kDigitizerBins, kDigitizerChannels, 0, kDigitizerChannels);

    auto pos_sig_mdt_raw = TCAHistogram<TH1D>("pos_sig_mdt_raw", "Positive Signal Module Time (Raw Data);ADC;Counts/Bin", kDigitizerBins, 0, kDigitizerBins);
    auto pos_sig_trt_raw = TCAHistogram<TH2D>("pos_sig_trt_raw", "Positive Signal Trigger Time (Raw Data);ADC;Trigger ID;Counts/Bin", kDigitizerBins, 0, kDigitizerBins, 2, 0, 2);

    // Calibrated Hists
    auto pos_sig_E = TCAHistogram<TH2D>("pos_sig_E", "Pos Sig Energy;Energy (keV),Crystal,Counts/Bin", kMaxEnergy / kEnergyPerBin, 0, kMaxEnergy, kDigitizerChannels, 0, kDigitizerChannels);
    auto pos_sig_cht = TCAHistogram<TH2D>("pos_sig_cht", "Pos Sig Channel Time;Channel Time (ns);Channel;Counts/Bin", kDigitizerBins, 0, (kDigitizerBins)*kNsPerBin, kDigitizerChannels, 0, kDigitizerChannels);
    auto pos_sig_mdt = TCAHistogram<TH1D>("pos_sig_mdt", "Positive Signal Module Time;Module Time (ns);Counts/Bin", kDigitizerBins, 0, (kDigitizerBins)*kNsPerBin);
    auto pos_sig_trt = TCAHistogram<TH2D>("pos_sig_trt", "Positive Signal Trigger Time;Trigger Time (ns);Trigger ID;Counts/Bin", kDigitizerBins, 0, (kDigitizerBins)*kNsPerBin, 2, 0, 2);

    // Sum Hists
    auto pos_sig_sum = TCAHistogram<TH2D>("pos_sig_sum", "Pos Sig Energy (Detector Sum);Energy (keV);Detector;Counts/bin", kMaxEnergy / kEnergyPerBin, 0, kMaxEnergy, kDigitizerChannels / 4, 0, kDigitizerChannels / 4);

    // Addback Hists
    auto pos_sig_addback = TCAHistogram<TH2D>("pos_sig_addback", "Pos Sig Energy (Detctor Addback);Energy (keV);Detetcor;Counts/Bin", kMaxEnergy / kEnergyPerBin, 0, kMaxEnergy, kDigitizerChannels / 4, 0, kDigitizerChannels / 4);
    auto pos_sig_addback_mult = TCAHistogram<TH1D>("pos_sig_addback_mult", "Pos Sig Addback Multiplicity;Multiplicity;Counts/Bin", 4, 0, 4);

    #endif // PROCESS_POS_SIG

    #if PROCESS_CEBR_ALL

    // Raw Hists
    auto cebr_all_inl_raw = TCAHistogram<TH2D>("cebr_all_inl_raw", "CeBr All Integration Long (Raw Data);ADC;Channel;Counts/Bin", kDigitizerBins, 0, kDigitizerBins, kDigitizerChannels, 0, kDigitizerChannels);
    auto cebr_all_cht_raw = TCAHistogram<TH2D>("cebr_all_cht_raw", "CeBr All Channel Time (Raw Data);ADC;Channel;Counts/Bin", kDigitizerBins, 0, kDigitizerBins, kDigitizerChannels, 0, kDigitizerChannels);
    auto cebr_all_ins_raw = TCAHistogram<TH2D>("cebr_all_ins_raw", "CeBr All Integration Short (Raw Data);ADC;Channel;Counts/Bin", kDigitizerBins, 0, kDigitizerBins, kDigitizerChannels, 0, kDigitizerChannels);

    auto cebr_all_mdt_raw = TCAHistogram<TH1D>("cebr_all_mdt_raw", "CeBr All Module Time (Raw Data);ADC;Counts/Bin", kDigitizerBins, 0, kDigitizerBins);
    auto cebr_all_trt_raw = TCAHistogram<TH2D>("cebr_all_trt_raw", "CeBr All Trigger Time (Raw Data);ADC;Trigger ID;Counts/Bin", kDigitizerBins, 0, kDigitizerBins, 2, 0, 2);

    // Calibrated Hists
    auto cebr_all_El = TCAHistogram<TH2D>("cebr_all_El", "CeBr All Energy (Long Integration);Energy (keV);Channel;Counts/Bin", kMaxEnergy / kEnergyPerBin, 0, kMaxEnergy, kDigitizerChannels, 0, kDigitizerChannels);
    auto cebr_all_cht = TCAHistogram<TH2D>("cebr_all_cht", "CeBr All Channel Time; Channel Time (ns);Channel;Counts/Bin", kDigitizerBins, 0, (kDigitizerBins)*kNsPerBin, kDigitizerChannels, 0, kDigitizerChannels);
    auto cebr_all_Es = TCAHistogram<TH2D>("cebr_all_Es", "CeBr All Energy (Short Integration); Energy (keV);Channel;Counts/Bin", kDigitizerBins, 0, kDigitizerBins, kDigitizerChannels, 0, kDigitizerChannels);

    auto cebr_all_mdt = TCAHistogram<TH1D>("cebr_all_mdt", "CeBr All Module Time;Module Time (ns);Counts/Bin", kDigitizerBins, 0, (kDigitizerBins)*kNsPerBin);
    auto cebr_all_trt = TCAHistogram<TH2D>("cebr_all_trt", "CeBr All Trigger Time;Trigger Time (ns);Trigger ID;Counts/Bin", kDigitizerBins, 0, (kDigitizerBins)*kNsPerBin, 2, 0, 2);

    #endif // PROCESS_CEBR_ALL

} // namespace CAHistograms

#endif // CAHISTOGRAMS_HPP