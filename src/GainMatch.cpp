// C++ includes
#include <iostream>
#include <memory>

// ROOT Includes
#include <TFile.h>
#include <TH2D.h>
#include <TSpectrum.h>

// Project Includes

namespace Configuration
{
    inline const char *kAmplitudeHistogramName = "clover_cross/cc_amp"; // Path to the amplitude histogram within the root file

    inline constexpr unsigned int kBackgroundSmoothing = 50; // Smoothing parameter for background subtraction
    inline const char *kBackgroundOptions = "Compton";       // Sigma for peak finding
    inline constexpr unsigned int kRebinFactor = 4;          // Rebin factor for histograms

    inline constexpr unsigned int kMaxPeaks = 10;                                                // Maximum number of peaks to identify
    inline constexpr std::pair<double, double> kPeakSearchRange = std::make_pair(6000, 22000.0); // Bin range in which to search for peaks
    inline constexpr double kPeakSigma = 15.0;                                                   // Expected peak width (sigma) in bins
    inline constexpr double kPeakThreshold = 0.10;                                               // Minimum height (as fraction of max) to consider a peak
    inline constexpr double kPeakCentroidRatio = 2614.511 / 1460.820;                            // Ratio describing 208Tl and 40K background lines
    inline constexpr double kPeakRatioTolerance = 0.003;                                         // Acceptable deviation in ratio when matching peaks

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

    std::cout << "Gain matching " << input_filename << " to " << reference_filename << " and saving to " << output_filename << std::endl;

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

    TSpectrum spectrum_utility(kMaxPeaks);

    // Perform background subtraction to give peak finding a better chance

    for (auto hist : {reference_hist, input_hist})
    {
        for (size_t i = 1; i <= hist->GetNbinsY(); ++i) // Loop over each bin in Y (channel)
        {
            auto *proj = hist->ProjectionX("_px", i, i);
            auto *bg = spectrum_utility.Background(proj, kBackgroundSmoothing, kBackgroundOptions);
            proj->Add(bg, -1);
            for (int j = 0; j <= proj->GetNbinsX(); ++j) // Copy back to original 2D histogram
                hist->SetBinContent(i, j, proj->GetBinContent(j));
            delete proj;
            delete bg;
        }
    }

    // Find Peaks
    // Recall that we rebinned the histograms, so the peak positions need to be scaled back
    std::vector<std::vector<double>> reference_peaks, input_peaks;
    for (auto &[hist, peaks] : {std::make_pair(reference_hist, std::ref(reference_peaks)),
                                std::make_pair(input_hist, std::ref(input_peaks))}) // Loop over both histograms
    {
        bool ref = (&hist == &reference_hist);
        for (size_t i = 1; i <= hist->GetNbinsY(); ++i) // Loop over each channel
        {
            auto hist_proj = hist->ProjectionX("_px", i, i);
            hist_proj->GetXaxis()->SetRange(kPeakSearchRange.first / kRebinFactor, kPeakSearchRange.second / kRebinFactor); // Set search range
            auto nfound = spectrum_utility.Search(hist_proj, kPeakSigma, "", kPeakThreshold);
            printf("%s: Found %u peaks in channel %zu (bin %zu):\n", ref ? "Reference" : "Input", nfound, i - 1, i); // Channel 0 is in bin 1
            auto xpeaks = spectrum_utility.GetPositionX();
            std::vector<double> channel_peaks;
            for (unsigned int p = 0; p < nfound; ++p)
            {
                printf("%f\n", xpeaks[p]);
                channel_peaks.push_back(xpeaks[p]);
            }
            peaks.push_back(channel_peaks);
            delete hist_proj;
        }
    }

    // Close input files
    reference_file->Close();
    input_file->Close();

    // Loop over channels and determine gain match parameters
    for (size_t ch = 0; ch < reference_peaks.size(); ++ch)
    {
        const auto &ref_ch_peaks = reference_peaks[ch];
        const auto &in_ch_peaks = input_peaks[ch];

        if (ref_ch_peaks.size() < 2 || in_ch_peaks.size() < 2)
        {
            std::cerr << "Not enough peaks found in channel " << ch << " to perform gain matching." << std::endl;
            continue;
        }

        // Find peaks that have the expected centroid ratio
        for (size_t i = 0; i < ref_ch_peaks.size(); ++i)
        {
            for (size_t j = 0; j < i; ++j)
            {
                double ref_ratio = ref_ch_peaks[i] / ref_ch_peaks[j];
                if (std::abs(ref_ratio - kPeakCentroidRatio) / kPeakCentroidRatio < kPeakRatioTolerance) // Allow 5% tolerance
                {
                    printf("Channel %zu: Found matching peaks: Ref (%f,%f) with ratio match (%f)\n", ch, ref_ch_peaks[i], ref_ch_peaks[j], ref_ratio / kPeakCentroidRatio);
                }
            }
        }

        // Calculate gain and offset
        // double gain = (ref_peak2 - ref_peak1) / (in_peak2 - in_peak1);
        // double offset = ref_peak1 - gain * in_peak1;

        // std::cout << "Channel " << ch << ": Gain = " << gain << ", Offset = " << offset << std::endl;
    }

    return 0;
}
