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
#include "Configuration.hpp"

/* #endregion Includes */

/* #region DAQ Channel Key */
// Define which DAQ channels correspond to which detectors (C = Clover HPGe, S = Single Crystal HPGe, c = CeBr, l = LaBr) see elog for details

// clover_cross
#define C1E1 0
#define C1E2 1
#define C1E3 2
#define C1E4 3
#define C3E1 4
#define C3E2 5
#define C3E3 6
#define C3E4 7
#define C5E1 8
#define C5E2 9
#define C5E3 10
#define C5E4 11
#define C7E1 12
#define C7E2 13
#define C7E3 14
#define C7E4 15

// clover_back
#define B1E1 0
#define B1E2 1
#define B1E3 2
#define B1E4 3
#define B2E1 4
#define B2E2 5
#define B2E3 6
#define B2E4 7
#define B3E1 8
#define B3E2 9
#define B3E3 10
#define B3E4 11
#define B5E1 12
#define B5E2 13
#define B5E3 14
#define B5E4 15

// pos_sig
#define ZDEG 0
#define SB4E1 2
#define B4E1 4
#define B4E2 5
#define B4E3 6
#define B4E4 7

// cebr_all
#define cB 0
#define cC 1
#define cD 2
#define cF 3
#define cG 4
#define cH 5
#define cK 6
#define cO 7
#define cBJ 8
#define cBK 9
#define cBL 10
#define L3 11
#define MPAD 12

static const std::map<std::string, int> cross_channel_map = {
    {"C1E1", C1E1}, {"C1E2", C1E2}, {"C1E3", C1E3}, {"C1E4", C1E4}, {"C3E1", C3E1}, {"C3E2", C3E2}, {"C3E3", C3E3}, {"C3E4", C3E4}, {"C5E1", C5E1}, {"C5E2", C5E2}, {"C5E3", C5E3}, {"C5E4", C5E4}, {"C7E1", C7E1}, {"C7E2", C7E2}, {"C7E3", C7E3}, {"C7E4", C7E4}};
static const std::map<std::string, int> back_channel_map = {
    {"B1E1", B1E1}, {"B1E2", B1E2}, {"B1E3", B1E3}, {"B1E4", B1E4}, {"B2E1", B2E1}, {"B2E2", B2E2}, {"B2E3", B2E3}, {"B2E4", B2E4}, {"B3E1", B3E1}, {"B3E2", B3E2}, {"B3E3", B3E3}, {"B3E4", B3E4}, {"B5E1", B5E1}, {"B5E2", B5E2}, {"B5E3", B5E3}, {"B5E4", B5E4}};
static const std::map<std::string, int> possig_channel_map = {
    {"ZDEG", ZDEG}, {"SB4E1", SB4E1}, {"B4E1", B4E1}, {"B4E2", B4E2}, {"B4E3", B4E3}, {"B4E4", B4E4}};
static const std::map<std::string, int> cebr_channel_map = {
    {"cB", cB}, {"cC", cC}, {"cD", cD}, {"cF", cF}, {"cG", cG}, {"cH", cH}, {"cK", cK}, {"cO", cO}, {"cBJ", cBJ}, {"cBK", cBK}, {"cBL", cBL}, {"L3", L3}, {"MPAD", MPAD}};

/* #endregion DAQ Channel Key*/

/* #region Helper Functions*/

// Add-back function
double cloverAddBackEnergy(std::vector<double> crystal_energies, std::vector<double> crystal_times)
{
    // Clover add-back function, modified from Samantha's code

    if (crystal_energies.size() != crystal_times.size())
    {
        std::cerr << "Error: cloverAddBackEnergy requires the same number of detector energies and times" << std::endl;
        return 0;
    }

    double highest_energy = 0;
    int highest_energy_index = -1;
    double final_energy = 0;
    std::vector<double> delta_t;
    unsigned short mult = 0; // Addback multiplicity counter

    // Find highest energy hit, this is most likely the primary hit
    for (size_t crystal = 0; crystal < crystal_energies.size(); crystal++)
    {
        if (crystal_energies[crystal] > CLOVER_ADDBACK_THRESHOLD) // Energy cut in keV
        {
            if (crystal_energies[crystal] > highest_energy)
            {
                highest_energy = crystal_energies[crystal];
                highest_energy_index = crystal;
            }
        }
    }
    if (highest_energy_index > -1) // If a primary hit was found
    {
        for (size_t crystal = 0; crystal < crystal_energies.size(); crystal++) // Perform the addback
        {
            delta_t.push_back(fabs(Histograms::kNsPerBin * (crystal_times[highest_energy_index] - crystal_times[crystal])));
            if (fabs(Histograms::kNsPerBin * (crystal_times[highest_energy_index] - crystal_times[crystal])) < CLOVER_ADDBACK_WINDOW) // Time cut in ns
            {
                if (crystal_energies[crystal] > CLOVER_ADDBACK_THRESHOLD) // Energy cut in keV
                {
                    final_energy += crystal_energies[crystal];
                    mult++;
                }
            }
        }
    }

    // fmt::print("Energies: {}\n Times: {}\n Final Energy: {}\n\n", crystal_energies, delta_t, final_energy);

    return final_energy;
}

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

TSpline3 createCalSpline(const std::string &filename)
{
    std::vector<double> knot_x, knot_y;
    std::ifstream calfile(filename);
    if (!calfile.is_open())
    {
        std::cerr << "Failed to open calibration file: " << filename << std::endl;
        return TSpline3();
    }

    // Read spline knots from the calibration file
    std::string line;
    while (std::getline(calfile, line))
    {
        bool skipped_linear = false;
        // Helper lambda to trim whitespace
        auto trim = [](std::string &s)
        {
            const char *whitespace = " \t\n\r\f\v";
            size_t start = s.find_first_not_of(whitespace);
            if (start == std::string::npos)
            {
                s.clear();
                return;
            }
            size_t end = s.find_last_not_of(whitespace);
            s = s.substr(start, end - start + 1);
        };

        while (std::getline(calfile, line))
        {
            trim(line);
            if (line.empty())
                continue;
            if (line[0] == '#')
                continue;

            // The first non-comment/non-empty line contains the linear parameters (slope offset).
            // Skip that line and only start collecting knots afterwards.
            if (!skipped_linear)
            {
                skipped_linear = true;
                continue;
            }

            std::istringstream iss(line);
            double x, y;
            if (iss >> x >> y)
            {
                knot_x.push_back(x);
                knot_y.push_back(y);
            }
        }
        if (line.empty() || line[0] == '#')
            continue;

        std::istringstream iss(line);
        double x, y;
        if (iss >> x >> y)
        {
            knot_x.push_back(x);
            knot_y.push_back(y);
        }
    }

    calfile.close();

    if (knot_x.empty())
    {
        std::cerr << "No spline knots found in " << filename << std::endl;
        return TSpline3();
    }

    // Create and return the spline
    return TSpline3("spline", knot_x.data(), knot_y.data(), knot_x.size(), "b1e1");
}

std::vector<double> readLinearCalParams(const std::string &filename)
{
    std::vector<Double_t> params(2, 0.0); // Initialize with two zeros
    std::ifstream calfile(filename);
    if (!calfile.is_open())
    {
        std::cerr << "Failed to open calibration file: " << filename << std::endl;
        return params;
    }

    std::string line;
    while (std::getline(calfile, line))
    {
        // Trim leading whitespace
        size_t first = line.find_first_not_of(" \t\r\n");
        if (first == std::string::npos)
            continue; // empty line

        // Skip comment lines
        if (line[first] == '#')
            continue;

        std::istringstream iss(line);
        double offset, slope;
        if (iss >> offset >> slope)
        {
            params[0] = offset;
            params[1] = slope;
            break;
        }
        // If the line didn't contain two valid numbers, keep searching
    }

    calfile.close();
    return params;
}

std::function<double(double)> makeCalibration(std::vector<double> linear_params, TSpline3 cal_spline)
{

    double offset = linear_params[0];
    double slope = linear_params[1];

    // Create calibration function: output = slope * input + offset + spline(input)
    auto calibration_func = [slope, offset, cal_spline](double input) -> double
    {
        double lincal_E = slope * input + offset;
        double spline_corr = cal_spline.Eval(lincal_E);
        double energy = lincal_E + spline_corr;
        // std::cout << "Input: " << input << ", lincal E: " << lincal_E << ", Spline Corr: " << spline_corr << ", Final E: " << energy << std::endl;
        return energy;
    };

    // Return the callable calibration lambda instead of an empty std::function
    return calibration_func;
}

/* #endregion Helper Functions*/

// Main function
int main(int argc, char *argv[])
{
    if (argc != 6)
    {
        std::cerr << "Usage: " << argv[0] << " <calibration directory> <gain shift directory> <run file directory> <run number> <output filename>" << std::endl;
        return 1;
    }

    std::cout << "Welcome to CASort!" << std::endl;

    std::string calibration_dir = argv[1];
    std::cout << "Using calibration directory: " << calibration_dir << std::endl;

    std::string gain_shift_dir = argv[2];
    std::cout << "Using gain shift directory: " << gain_shift_dir << std::endl;

    std::string run_file_dir = argv[3];
    std::string run_number = Form("%03d", std::stoi(argv[4]));
    std::string run_file_name = Form(RUN_FILE_NAME_TEMPLATE, run_number.c_str());
    std::string input_filename = Form("%s/%s", run_file_dir.c_str(), run_file_name.c_str());
    std::cout << "Input file: " << input_filename << std::endl;

    std::string output_filename = argv[5];

    // Set number of threads
    std::cout << "Using " << N_THREADS << " threads!" << std::endl;
    ROOT::EnableImplicitMT(N_THREADS);
    ROOT::EnableThreadSafety();

    /* #region Calibration Setup */

    // Open calibration file and make calibration functions
    std::vector<std::function<double(double)>> cloverCrossECal, cloverBackECal, posSigECal, cebrAllECal;

    // Load calibration splines and linear params

    // Clover Cross
#if PROCESS_CLOVER_CROSS
    std::vector<TSpline3> clover_cross_cal_splines;
    std::vector<std::vector<Double_t>> clover_cross_cal_linear_params;

    for (int det : {1, 3, 5, 7})
    {
        for (int xtal = 1; xtal <= 4; xtal++)
        {
            std::string cal_filename = Form("%s/C%iE%i.cal_params.txt", calibration_dir.c_str(), det, xtal);
            cloverCrossECal.push_back(makeCalibration(readLinearCalParams(cal_filename), createCalSpline(cal_filename)));
        }
    }
#endif // PROCESS_CLOVER_CROSS

    // Clover Back
#if PROCESS_CLOVER_BACK
    std::vector<TSpline3> clover_back_cal_splines;
    std::vector<std::vector<Double_t>> clover_back_cal_linear_params;
    for (int det : {1, 2, 3, 5})
    {
        for (int xtal = 1; xtal <= 4; xtal++)
        {
            std::string cal_filename = Form("%s/B%iE%i.cal_params.txt", calibration_dir.c_str(), det, xtal);
            cloverBackECal.push_back(makeCalibration(readLinearCalParams(cal_filename), createCalSpline(cal_filename)));
        }
    }
#endif // PROCESS_CLOVER_BACK

    /* #endregion Calibration Setup */

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
    auto fillHistograms = [&](TTreeReader &eventReader)
    {
        /* #region Set the branch addresses for the TTree */

        TTreeReaderArray<double> cc_amp_val(eventReader, "clover_cross.amplitude");
        TTreeReaderArray<double> cc_cht_val(eventReader, "clover_cross.channel_time");
        TTreeReaderArray<double> cc_mdt_val(eventReader, "clover_cross.module_timestamp");
        TTreeReaderArray<double> cc_plu_val(eventReader, "clover_cross.pileup");
        TTreeReaderArray<double> cc_trt_val(eventReader, "clover_cross.trigger_time");

        TTreeReaderArray<double> cb_amp_val(eventReader, "clover_back.amplitude");
        TTreeReaderArray<double> cb_cht_val(eventReader, "clover_back.channel_time");
        TTreeReaderArray<double> cb_mdt_val(eventReader, "clover_back.module_timestamp");
        TTreeReaderArray<double> cb_plu_val(eventReader, "clover_back.pileup");
        TTreeReaderArray<double> cb_trt_val(eventReader, "clover_back.trigger_time");

        TTreeReaderArray<double> ps_amp_val(eventReader, "pos_sig.amplitude");
        TTreeReaderArray<double> ps_cht_val(eventReader, "pos_sig.channel_time");
        TTreeReaderArray<double> ps_mdt_val(eventReader, "pos_sig.module_timestamp");
        TTreeReaderArray<double> ps_plu_val(eventReader, "pos_sig.pileup");
        TTreeReaderArray<double> ps_trt_val(eventReader, "pos_sig.trigger_time");

        TTreeReaderArray<double> ce_inl_val(eventReader, "cebr_all.integration_long");
        TTreeReaderArray<double> ce_cht_val(eventReader, "cebr_all.channel_time");
        TTreeReaderArray<double> ce_mdt_val(eventReader, "cebr_all.module_timestamp");
        TTreeReaderArray<double> ce_ins_val(eventReader, "cebr_all.integration_short");
        TTreeReaderArray<double> ce_trt_val(eventReader, "cebr_all.trigger_time");

        /* #endregion */

        /* #region Get Histogram pointers*/

        // Use histograms defined in Histograms.hpp

        // Clover Cross Histograms

#if PROCESS_CLOVER_CROSS
        auto cc_amp = Histograms::cc_amp.GetThreadLocalPtr();
        auto cc_cht = Histograms::cc_cht.GetThreadLocalPtr();
        auto cc_plu = Histograms::cc_plu.GetThreadLocalPtr();
        auto cc_trt = Histograms::cc_trt.GetThreadLocalPtr();
        auto cc_mdt = Histograms::cc_mdt.GetThreadLocalPtr();
        auto cc_E = Histograms::cc_E.GetThreadLocalPtr();
        auto cc_sum = Histograms::cc_sum.GetThreadLocalPtr();
        auto cc_abE = Histograms::cc_abE.GetThreadLocalPtr();
        auto cc_abM = Histograms::cc_abM.GetThreadLocalPtr();
#endif // PROCESS_CLOVER_CROSS

#if PROCESS_CLOVER_BACK

#endif // PROCESS_CLOVER_BACK

#if PROCESS_POS_SIG

#endif // PROCESS_POS_SIG

#if PROCESS_CEBR_ALL

#endif // PROCESS_CEBR_ALL

        /* #endregion */

        // Loop over the entries in the tree
        while (eventReader.Next())
        {

#if PROCESS_CLOVER_CROSS

            // Module Time
            cc_mdt->Fill(cc_mdt_val[0] * Histograms::kNsPerBin);

            // Trigger Times
            cc_trt->Fill(cc_trt_val[0] * Histograms::kNsPerBin, 0);
            cc_trt->Fill(cc_trt_val[1] * Histograms::kNsPerBin, 1);
            // Main Loop

            // Detector Loop
            for (size_t det = 0; det < 4; det++)
            {
                std::vector<double> xtal_E;
                std::vector<double> xtal_T;

                // Crystal Loop
                for (size_t xtal = 0; xtal < 4; xtal++)
                {
                    auto ch = det * 4 + xtal; // Channel number 0-15

                    // Raw Histograms
                    cc_amp->Fill(cc_amp_val[ch], ch);
                    cc_cht->Fill(cc_cht_val[ch], ch);
                    cc_plu->Fill(cc_plu_val[ch], ch);

                    // Calibrated Histograms
                    if (!std::isnan(cc_amp_val[ch]) && !std::isnan(cc_cht_val[ch]))
                    {
                        // std::cout << "Channel: " << ch << ", ";
                        double energy = cloverCrossECal[ch](cc_amp_val[ch]);
                        double cht = cc_cht_val[ch] * Histograms::kNsPerBin;
                        cc_E->Fill(energy, ch);
                        cc_cht->Fill(cht, ch);
                        cc_sum->Fill(energy, det); // ch / 4 is the detector number
                        xtal_E.push_back(energy);
                        xtal_T.push_back(cht);
                    }
                }

                if (!xtal_E.empty())
                {
                    cc_abE->Fill(cloverAddBackEnergy(xtal_E, xtal_T), det);
                    cc_abM->Fill(xtal_E.size(), det);
                }
            }

#endif // PROCESS_CLOVER_CROSS

#if PROCESS_CLOVER_BACK

#endif // PROCESS_CLOVER_BACK

#if PROCESS_POS_SIG

#endif // PROCESS_POS_SIG

#if PROCESS_CEBR_ALL

#endif // PROCESS_CEBR_ALL

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

// Clover Cross Histograms
#if PROCESS_CLOVER_CROSS
    auto cc_dir = outfile->mkdir("clover_cross");
    cc_dir->cd();
    Histograms::cc_amp.Write();
    Histograms::cc_cht.Write();
    Histograms::cc_plu.Write();
    Histograms::cc_trt.Write();
    Histograms::cc_mdt.Write();
    Histograms::cc_E.Write();
    Histograms::cc_sum.Write();
    Histograms::cc_abE.Write();
    Histograms::cc_abM.Write();

    outfile->cd();
#endif // PROCESS_CLOVER_CROSS

// Clover Back Histograms
#if PROCESS_CLOVER_BACK

#endif // PROCESS_CLOVER_BACK

// Positive Signal Histograms
#if PROCESS_POS_SIG

#endif // PROCESS_POS_SIG

// CeBr All Histograms
#if PROCESS_CEBR_ALL

#endif // PROCESS_CEBR_ALL

    /* #endregion */

    std::cout << "Saved histograms to file: " << outfile->GetName() << std::endl;

    outfile->Close();
    delete outfile;

    std::cout << "Done!" << std::endl;

    return 0;
}