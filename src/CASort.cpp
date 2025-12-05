/**
 * @file extract_hists.cpp
 * @brief This program extracts histograms from a ROOT file and saves them to a new ROOT file.
 *
 * This program takes a ROOT file as input, extracts a TTree from it, and creates two histograms
 * from the data in the TTree. The histograms are then saved to a new ROOT file.
 *
 * Compile with: g++ `root-config --cflags --libs` -o extract_hists extract_hists.cpp
 * (May need to put flags after file arguments on Ubuntu)
 *
 * Usage:
 * @code
 * ./extract_hists <filename>
 * @endcode
 *
 * @param argc The number of command-line arguments.
 * @param argv The array of command-line arguments.
 * @return Returns 0 on success, or 1 on failure.
 */

/* #region Includes */

// C++ Includes
#include <vector>
#include <string>
#include <fstream>
#include <sstream>
#include <atomic>
#include <thread>
#include <chrono>
#include <filesystem>

// ROOT Includes
#include <TFile.h>
#include <TTree.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TStopwatch.h>
#include <TTreeReader.h>
#include <TTreeReaderArray.h>
#include <ROOT/TThreadedObject.hxx>
#include <ROOT/TTreeProcessorMT.hxx>
#include <TSpline.h>

// Project Includes
#include "Histograms.hpp"

/* #endregion Includes */

/* #region Global Config */
// Number of Hardware Threads
#define N_THREADS 30 // Number of threads to use for processing, defaults to system max

/* #endregion Global Config */

/* #region Helper Functions*/

void displayProgressBar(std::atomic<ULong64_t> &processedEntries, ULong64_t totalEntries)
{
    const int barWidth = 50; // Width of the progress bar
    while (processedEntries < totalEntries)
    {
        double progress = static_cast<double>(processedEntries) / totalEntries;
        int pos = static_cast<int>(barWidth * progress);

        std::cout << "[";
        for (int i = 0; i < barWidth; ++i)
        {
            if (i < pos)
                std::cout << "=";
            else if (i == pos)
                std::cout << ">";
            else
                std::cout << " ";
        }
        std::cout << "] " << int(progress * 100.0) << "% (" << processedEntries << "/" << totalEntries << ")\r";
        std::cout.flush();

        std::this_thread::sleep_for(std::chrono::milliseconds(100)); // Update every 100ms
    }
    std::cout << "[";
    for (int i = 0; i < barWidth; ++i)
        std::cout << "=";
    std::cout << "] 100% (" << totalEntries << "/" << totalEntries << ")\n";
}

/* #endregion Helper Functions*/

// Main function
int main(int argc, char *argv[])
{
    if (argc != 3)
    {
        std::cerr << "Usage: " << argv[0] << " <input filename> <output filename>" << std::endl;
        return 1;
    }
    std::string input_filename = argv[1];
    std::string output_filename = argv[2];

    // Set number of threads
    std::cout << "Started extract_hists using " << N_THREADS << " threads!" << std::endl;
    ROOT::EnableImplicitMT(N_THREADS);

    /* #region Event Loop Setup*/

    auto infile = TFile::Open(input_filename.c_str());
    if (!infile || infile->IsZombie())
    {
        std::cerr << "Error opening input file" << std::endl;
        return 1;
    }
    std::cout << "Opened file: " << input_filename << std::endl;

    // Peak at the TTree to get the number of events
    TTree *tree;
    infile->GetObject("clover", tree);
    if (!tree)
    {
        std::cerr << "Error opening TTree" << std::endl;
        return 1;
    }

    ULong64_t n_entries = tree->GetEntries();

    std::cout << "Opened TTree: clover and counted " << n_entries << " events" << std::endl;

    delete tree;
    infile->Close();

    // Atomic counter for processed entries
    std::atomic<ULong64_t> processedEntries(0);

    // Start the progress bar in a separate thread
    std::thread progressBarThread(displayProgressBar, std::ref(processedEntries), n_entries);

    /* #endregion Event Loop Setup */

    std::cout << "Processing events..." << std::endl;

    // Create a TTreeReader to read the TTree
    ROOT::TTreeProcessorMT EventProcessor(input_filename.c_str(), "clover");

    // Fill Function
    using namespace Histograms;
    auto fillHistograms = [&](TTreeReader &eventReader)
    {
        /* #region Set the branch addresses for the TTree */

        TTreeReaderArray<double> cc_amp(eventReader, "clover_cross.amplitude");
        TTreeReaderArray<double> cc_cht(eventReader, "clover_cross.channel_time");
        TTreeReaderArray<double> cc_mdt(eventReader, "clover_cross.module_timestamp");
        TTreeReaderArray<double> cc_plu(eventReader, "clover_cross.pileup");
        TTreeReaderArray<double> cc_trt(eventReader, "clover_cross.trigger_time");

        TTreeReaderArray<double> cb_amp(eventReader, "clover_back.amplitude");
        TTreeReaderArray<double> cb_cht(eventReader, "clover_back.channel_time");
        TTreeReaderArray<double> cb_mdt(eventReader, "clover_back.module_timestamp");
        TTreeReaderArray<double> cb_plu(eventReader, "clover_back.pileup");
        TTreeReaderArray<double> cb_trt(eventReader, "clover_back.trigger_time");

        TTreeReaderArray<double> ps_amp(eventReader, "pos_sig.amplitude");
        TTreeReaderArray<double> ps_cht(eventReader, "pos_sig.channel_time");
        TTreeReaderArray<double> ps_mdt(eventReader, "pos_sig.module_timestamp");
        TTreeReaderArray<double> ps_plu(eventReader, "pos_sig.pileup");
        TTreeReaderArray<double> ps_trt(eventReader, "pos_sig.trigger_time");

        TTreeReaderArray<double> ce_inl(eventReader, "cebr_all.integration_long");
        TTreeReaderArray<double> ce_cht(eventReader, "cebr_all.channel_time");
        TTreeReaderArray<double> ce_mdt(eventReader, "cebr_all.module_timestamp");
        TTreeReaderArray<double> ce_ins(eventReader, "cebr_all.integration_short");
        TTreeReaderArray<double> ce_trt(eventReader, "cebr_all.trigger_time");

        /* #endregion */

        /* #region Get Histogram pointers*/

        // Use histograms defined in Histograms.hpp

        // Clover Cross Histograms
        auto cc_amp_raw_hist = clover_cross_amp_raw.MakePtr();
        auto cc_cht_raw_hist = clover_cross_cht_raw.MakePtr();
        auto cc_plu_hist = clover_cross_plu.MakePtr();
        auto cc_trt_raw_hist = clover_cross_trt_raw.MakePtr();
        auto cc_mdt_raw_hist = clover_cross_mdt_raw.MakePtr();

        auto cc_E_hist = clover_cross_E.MakePtr();
        auto cc_cht_hist = clover_cross_cht.MakePtr();
        auto cc_mdt_hist = clover_cross_mdt.MakePtr();
        auto cc_trt_hist = clover_cross_trt.MakePtr();

        auto cc_sum_hist = clover_cross_sum.MakePtr();
        auto cc_adb_hist = clover_cross_addback.MakePtr();
        auto cc_adb_mult_hist = clover_cross_addback_mult.MakePtr();

        /* #endregion */

        // Loop over the entries in the tree
        while (eventReader.Next())
        {
            /* #region clover_cross */

            // Module Time
            cc_mdt_raw_hist->Fill(cc_mdt[0]);
            cc_mdt_hist->Fill(cc_mdt[0] * kNsPerBin);

            // Trigger Times
            cc_trt_raw_hist->Fill(cc_trt[0], 0);
            cc_trt_hist->Fill(cc_trt[0] * kNsPerBin, 0);
            cc_trt_raw_hist->Fill(cc_trt[1], 1);
            cc_trt_hist->Fill(cc_trt[1] * kNsPerBin, 1);

            // Main Loop

            // Detector Loop
            for (size_t det = 0; det < 4; det++)
            {
                std::vector<double> xtal_energies;
                std::vector<double> xtal_times;

                // Crystal Loop
                for (size_t xtal = 0; xtal < 4; xtal++)
                {
                    auto ch = det * 4 + xtal; // Channel number 0-15

                    // Raw Histograms
                    cc_amp_raw_hist->Fill(cc_amp[ch], ch);
                    cc_cht_raw_hist->Fill(cc_cht[ch], ch);
                    cc_plu_hist->Fill(cc_plu[ch], ch);

                    // Calibrated Histograms
                    if (!std::isnan(cc_amp[ch]) && !std::isnan(cc_cht[ch]))
                    {

                        double energy = cc_amp[ch] / 0.17; // cloverCrossECal[ch](cc_amp[ch]);
                        double cht = cc_cht[ch] * kNsPerBin;

                        cc_E_hist->Fill(energy, ch);
                        cc_cht_hist->Fill(cht, ch);
                        cc_sum_hist->Fill(energy, det); // ch / 4 is the detector number
                        xtal_energies.push_back(energy);
                        xtal_times.push_back(cht);
                    }
                }

                if (!xtal_energies.empty())
                {
                    cc_adb_hist->Fill(xtal_energies[0], det);
                    cc_adb_mult_hist->Fill(xtal_energies.size(), det);
                }
            }

            /* #endregion */

            processedEntries++;
        }
    };

    // Loop over the entries in the TTree and fill the histograms appropriately
    TStopwatch timer;
    timer.Start();
    EventProcessor.Process(fillHistograms);
    timer.Stop();

    progressBarThread.join();

    std::cout << "Processed events in " << timer.RealTime() << " seconds ("
              << static_cast<double>(processedEntries) / timer.RealTime()
              << " events/second)" << std::endl;

    // Save the histograms to a new ROOT file
    TFile *outfile = new TFile(output_filename.c_str(), "RECREATE");
    if (!outfile || outfile->IsZombie())
    {
        std::cerr << "Error creating output file" << std::endl;
        return 1;
    }

    /* #region Write Histograms */

    // Use histograms defined in Histograms.hpp
    using namespace Histograms;

    // Clover Cross Histograms
    auto cc_dir = outfile->mkdir("clover_cross");
    cc_dir->cd();
    clover_cross_amp_raw.Write();
    outfile->cd();

    /* #endregion */

    std::cout << "Saved histograms to file: " << outfile->GetName() << std::endl;

    outfile->Close();
    delete outfile;

    std::cout << "Done!" << std::endl;

    return 0;
}