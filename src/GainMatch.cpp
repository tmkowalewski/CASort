// C++ includes
#include <iostream>
#include <memory>

// ROOT Includes
#include <TFile.h>
#include <TH2D.h>
#include <TSpectrum.h>
#include <TMath.h>
#include <TF1.h>
#include <TFitResult.h>
#include <TFitResultPtr.h>
#include "GainMatch.h"

// Project Includes

// Configuration
#define DEBUG 1

namespace Configuration
{
    // Histogram Config
    inline const char *kAmplitudeHistogramName = "clover_cross/cc_amp"; // Path to the amplitude histogram within the root file
    inline constexpr unsigned int kNumChannels = 16;                    // Number of channels in the histogram

    // Peak Finding Config
    inline constexpr unsigned int kBackgroundSmoothing = 50;                                            // Smoothing parameter for background subtraction
    inline const char *kBackgroundOptions = "Compton";                                                  // Sigma for peak finding
    inline constexpr unsigned int kRebinFactor = 4;                                                     // Rebin factor for histograms
    inline constexpr std::pair<double, double> kReferenceEnergies = std::make_pair(1460.820, 2614.511); // Number of iterations for background subtraction
    inline constexpr unsigned int kMaxPeaks = 10;                                                       // Maximum number of peaks to identify
    inline constexpr std::pair<double, double> kPeakSearchRange = std::make_pair(6000, 22000.0);        // Bin range in which to search for peaks
    inline constexpr double kPeakSigma = 15.0;                                                          // Expected peak width (sigma) in bins
    inline constexpr double kPeakThreshold = 0.10;                                                      // Minimum height (as fraction of max) to consider a peak
    inline constexpr double kPeakCentroidRatio = 2614.511 / 1460.820;                                   // Ratio describing 208Tl and 40K background lines
    inline constexpr double kPeakRatioTolerance = 0.0025;                                               // Acceptable deviation in ratio when matching peaks

    // Fitting Config
    inline constexpr double kFitBounds = 3;                        // How many sigma to include in fit range
    inline constexpr double kFitBgA = 0, kFitBgB = 0, kFitBgC = 0; // Quadratic background initial params Cx^2 + Bx + A
    inline constexpr double kFitR = 50;                            // Initial skew amplitude parameter R (%)
    inline constexpr double kFitBeta = 0.3;                        // Initial skew degree parameter beta
    inline constexpr double kFitBgStep = 0;                        // Initial background step parameter (%)
    inline constexpr double kFitFWHM = 15;                         // Initial FWHM

}

void BackgroundSubtraction(TH2D *reference_hist, TH2D *input_hist, TSpectrum *spectrum_utility)
{
    for (auto hist : {reference_hist, input_hist})
    {
        for (size_t ch = 0; ch < hist->GetNbinsY(); ++ch) // Loop over each bin in Y (channel)
        {
            auto *proj = hist->ProjectionX("_px", ch + 1, ch + 1);
            auto *bg = spectrum_utility->Background(proj, Configuration::kBackgroundSmoothing, Configuration::kBackgroundOptions);
            proj->Add(bg, -1);
            for (int j = 0; j <= proj->GetNbinsX(); ++j) // Copy back to original 2D histogram
                hist->SetBinContent(ch + 1, j, proj->GetBinContent(j));
            delete proj;
            delete bg;
        }
    }
}

int main(int argc, char *argv[])
{
    // Introduction
    if (argc != 4)
    {
        std::cerr << "Usage: " << argv[0] << " <reference file> <input file> <output file>" << std::endl;
        return 1;
    }

    const std::string reference_filename = argv[1];
    const std::string input_filename = argv[2];
    const std::string output_filename = argv[3];
    printf("============= Welcome to GainMatch! ===============\n");
    printf("------------- Current Configuration ---------------\n");
    printf("Using reference file: %s\n", reference_filename.c_str());
    printf("Using input file: %s\n", input_filename.c_str());
    printf("Output file: %s\n", output_filename.c_str());
    printf("---------------------------------------------------\n");

    std::cout << "Gain matching started!" << std::endl;

    using namespace Configuration;

    // Open Files

    auto reference_file = TFile::Open(reference_filename.c_str(), "READ");
    if (!reference_file)
    {
        std::cerr << "Error opening reference file" << std::endl;
        return 1;
    }

    auto input_file = TFile::Open(input_filename.c_str(), "READ");
    if (!input_file)
    {
        std::cerr << "Error opening input file" << std::endl;
        return 1;
    }

    // Get Histograms

    auto reference_hist = dynamic_cast<TH2D *>(reference_file->Get(kAmplitudeHistogramName));
    if (!reference_hist)
    {
        std::cerr << "Error retrieving reference histogram" << std::endl;
        return 1;
    }
    reference_hist->RebinX(kRebinFactor); // Rebin to make peakfinding easier

    auto input_hist = dynamic_cast<TH2D *>(input_file->Get(kAmplitudeHistogramName));
    if (!input_hist)

    {
        std::cerr << "Error retrieving input histogram" << std::endl;
        return 1;
    }
    input_hist->RebinX(kRebinFactor); // Rebin to make peakfinding easier

    // Make TSpectrum Object

    auto spectrum_utility = new TSpectrum(kMaxPeaks);

    // Perform background subtraction to give peak finding a better chance

    BackgroundSubtraction(reference_hist, input_hist, spectrum_utility);

    // Find Peaks and keep the ones that match the expected centroid ratio
    std::array<std::vector<double>, kNumChannels> reference_peaks, input_peaks;
    std::array<std::pair<double, double>, kNumChannels> matched_peaks_ref, matched_peaks_input; // Pair of peaks that look like reference energies

    for (auto &[hist, peaks] : {std::make_pair(reference_hist, std::ref(reference_peaks)),
                                std::make_pair(input_hist, std::ref(input_peaks))}) // Loop over both histograms
    {
        bool ref = (hist == reference_hist);

        // For each channel, find peaks
        for (size_t ch = 0; ch < hist->GetNbinsY(); ++ch) // Loop over each channel
        {
            auto hist_proj = hist->ProjectionX("_px", ch + 1, ch + 1);
            hist_proj->GetXaxis()->SetRange(kPeakSearchRange.first / kRebinFactor, kPeakSearchRange.second / kRebinFactor); // Set search range
            auto n_found = spectrum_utility->Search(hist_proj, kPeakSigma, "", kPeakThreshold);
            auto peak_pos = spectrum_utility->GetPositionX();
#if DEBUG >= 2
            printf("%s: Found %u peaks in channel %zu (bin %zu):\n", ref ? "Reference" : "Input", n_found, ch, ch + 1); // Channel 0 is in bin 1
#endif
            std::vector<double> channel_peaks;
            for (size_t p = 0; p < n_found; ++p)
            {
#if DEBUG >= 2
                printf("%.0f\n", peak_pos[p]); // For whatever reason we don't have to multiply by kRebinFactor here
#endif
                channel_peaks.push_back(peak_pos[p]);
            }
            std::sort(channel_peaks.begin(), channel_peaks.end());
            peaks[ch] = channel_peaks;
            delete hist_proj;
        }

        // Now find matching peaks based on expected centroid ratio
        for (size_t ch = 0; ch < peaks.size(); ++ch)
        {

            double best_ratio = 0.0;
            auto best_pair = std::make_pair(-1.0, -1.0);

            // Look through all combinations of peaks for this channel
            for (size_t i = 0; i < peaks[ch].size(); ++i)
            {
                for (size_t j = 0; j < i; ++j)
                {
                    auto ref_ratio = peaks[ch][i] / peaks[ch][j];
                    if (std::abs(ref_ratio - kPeakCentroidRatio) / kPeakCentroidRatio < kPeakRatioTolerance) // Allow tolerance about the expected ratio
                    {
                        if (std::abs(ref_ratio - kPeakCentroidRatio) < std::abs(best_ratio - kPeakCentroidRatio))
                        {
                            best_ratio = ref_ratio;
#if DEBUG >= 2
                            printf("%s: Found matching peaks: (%.0f,%.0f) in channel %zu with ratio match (%f)\n", ref ? "Reference" : "Input", peaks[ch][i], peaks[ch][j], ch, ref_ratio / kPeakCentroidRatio);
#endif
                        }

                        best_pair = std::make_pair(peaks[ch][j], peaks[ch][i]); // Store as (low, high)
                    }
                }
            }
            if (best_pair.first > 0 && best_pair.second > 0)
            {
                if (ref)
                    matched_peaks_ref[ch] = best_pair;
                else
                    matched_peaks_input[ch] = best_pair;
            }
            else
            {
                printf("%s: No matching peaks found within tolerance (%f) in channel %zu\n", ref ? "Reference" : "Input", kPeakRatioTolerance, ch);
            }
        }
    }

    // Actually fit the peaks at the locations we found to get more precise centroids
    std::array<std::pair<double, double>, kNumChannels> fitted_peaks_ref, fitted_peaks_input;
    for (size_t ch = 0; ch < kNumChannels; ++ch)
    {
        for (auto &[hist, matched_peaks] : {std::make_pair(reference_hist, std::ref(matched_peaks_ref)),
                                            std::make_pair(input_hist, std::ref(matched_peaks_input))})
        {
            bool ref = (hist == reference_hist);
            if (matched_peaks[ch].first < 0 || matched_peaks[ch].second < 0)
                continue; // No peaks found for this channel

            bool got_first = false;
            for (auto &peak_centroid : {std::ref(matched_peaks[ch].first), std::ref(matched_peaks[ch].second)})
            {
                auto hist_proj = hist->ProjectionX("_px", ch + 1, ch + 1);
                double peak_pos = peak_centroid.get();
                double peak_height = hist_proj->GetBinContent(hist_proj->FindBin(peak_pos));
                double fit_min = peak_pos - kFitBounds * kPeakSigma;
                double fit_max = peak_pos + kFitBounds * kPeakSigma;

                auto gaus = new TF1("gaus", "gaus", fit_min, fit_max);

                // Initial Values for Fit
                gaus->SetParameters(peak_height, peak_pos, kPeakSigma);
                gaus->SetParLimits(1, peak_pos - kPeakSigma, peak_pos + kPeakSigma);
                gaus->SetParLimits(2, 0, kPeakSigma * 2.0);

                hist_proj->Fit(gaus, "QLMRES0", "", fit_min, fit_max);
                auto mean = gaus->GetParameter(1);
#if DEBUG >= 2
                printf("%s: Fitted peak in channel %zu at initial pos %.1f to centroid %.3f\n", (hist == reference_hist) ? "Reference" : "Input", ch, peak_pos, mean);
#endif
                if (std::fabs(mean - peak_pos) > 5)
                {
                    std::cerr << (ref ? "Reference" : "Input") << ": Fit centroid (" << mean << ") deviated significantly from initial estimate (" << peak_pos << ") in channel " << ch << std::endl;
                }

                if (got_first == false)
                {
                    if (ref)
                        fitted_peaks_ref[ch].first = mean;
                    else
                        fitted_peaks_input[ch].first = mean;
                    got_first = true;
                }
                else
                {
                    if (ref)
                        fitted_peaks_ref[ch].second = mean;
                    else
                        fitted_peaks_input[ch].second = mean;
                }

                delete hist_proj;
                delete gaus;
            }
        }
    }

    // Loop over channels and determine gain match parameters
    for (size_t ch = 0; ch < kNumChannels; ++ch)
    {
        if (fitted_peaks_ref[ch].first < 0 || fitted_peaks_input[ch].first < 0)
        {
            std::cerr << "Channel " << ch << ": Missing fitted peaks, skipping gain match calculation" << std::endl;
            continue;
        }

        double ref_low = fitted_peaks_ref[ch].first;
        double ref_high = fitted_peaks_ref[ch].second;
        double input_low = fitted_peaks_input[ch].first;
        double input_high = fitted_peaks_input[ch].second;

#if DEBUG >= 1
        printf("Channel %zu: Fitted Reference Peaks: (%.3f, %.3f), Fitted Input Peaks: (%.3f, %.3f)\n", ch,
               fitted_peaks_ref[ch].first, fitted_peaks_ref[ch].second,
               fitted_peaks_input[ch].first, fitted_peaks_input[ch].second);
#endif

        // Calculate gain and offset
        double gain = (ref_high - ref_low) / (input_high - input_low);
        double offset = ref_low - gain * input_low;

        printf("Channel %zu: Gain = %.6f, Offset = %.10f\n", ch, gain, offset);
    }

    // Close input files
    reference_file->Close();
    delete reference_file;
    input_file->Close();
    delete input_file;

    return 0;
}
