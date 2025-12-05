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
#define N_THREADS std::thread::hardware_concurrency() // Number of threads to use for processing, defaults to system max

// Config Constants
#define CLOVER_ENERGY_THRESHOLD 150  // Energy threshold (keV) for clover detectors
#define CLOVER_ADDBACK_THRESHOLD 150 // Energy threshold (keV) for add-back
#define CLOVER_ADDBACK_WINDOW 150    // Time window (ns) around primary hit for add-back
#define ENERGY_PER_BIN 0.25          // Energy per bin (keV) for clover detectors for 16-bit digitizer
#define MAX_ENERGY 10000             // Maximum energy (keV) for clover detectors
#define CLOVER_COIN_WINDOW 100       // Time window (ns) for clover coincidence

#define POS_SIG_ENERGY_THRESHOLD 150 // Energy threshold (keV) for positive signal detectors

#define CEBR_ALL_ENERGY_THRESHOLD 150 // Energy threshold (keV) for CeBr (also LaBr and MPAD) detectors

// Calibration File Directories
#define ENERGY_CAL_DIR "energy_calibration/calibrations/"
#define GAIN_MATCH_DIR "energy_calibration/gain_matching/"

/* #endregion Global Config */

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

// Returns nominal energy given a vector of coefficients for a binomial of the form c[0]+c[1]*x
double calibrateLinear(double channel, std::vector<double> c)
{
    // Linear calibration function
    return c[0] + c[1] * channel;
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

std::vector<std::function<double(double)>> makeCalibrations(std::string calibrationDir)
{
    std::vector<std::function<double(double)>> calibrations;
    namespace fs = std::filesystem;

    try
    {
        // Check if directory exists
        if (!fs::exists(calibrationDir) || !fs::is_directory(calibrationDir))
        {
            std::cerr << "Calibration directory does not exist: " << calibrationDir << std::endl;
            return calibrations;
        }

        // Iterate through all files in the directory
        for (const auto &entry : fs::directory_iterator(calibrationDir))
        {
            if (!entry.is_regular_file())
                continue;

            std::string filename = entry.path().filename().string();

            // Check if file matches pattern "*.cal_params.txt"
            size_t pos = filename.find(".cal_params.txt");
            if (pos == std::string::npos || pos + 15 != filename.length())
                continue;

            std::string filepath = entry.path().string();
            std::ifstream calfile(filepath);

            if (!calfile.is_open())
            {
                std::cerr << "Failed to open calibration file: " << filepath << std::endl;
                continue;
            }

            // Read linear fit parameters (slope and offset)
            std::string line;
            double slope, offset;

            // Skip the comment line
            std::getline(calfile, line);

            // Read slope and offset
            if (!(calfile >> slope >> offset))
            {
                std::cerr << "Failed to read linear parameters from " << filepath << std::endl;
                calfile.close();
                continue;
            }

            // Skip the comment line about spline knots
            std::getline(calfile, line); // consume rest of previous line
            std::getline(calfile, line); // skip "# Spline knots: N"

            // Parse number of knots
            std::vector<double> knot_x, knot_y;
            std::string knot_line;

            while (std::getline(calfile, knot_line))
            {
                if (knot_line.empty() || knot_line[0] == '#')
                    continue;

                std::istringstream iss(knot_line);
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
                std::cerr << "No spline knots found in " << filepath << std::endl;
                continue;
            }

            // Create TSpline for this calibration file
            // We need to capture the spline in a way that persists with the lambda
            auto spline = std::make_shared<TSpline3>("spline",
                                                     knot_x.data(), knot_y.data(),
                                                     knot_x.size(), "b1e1");

            // Create calibration function: output = slope * input + offset + spline(input)
            auto calibration_func = [slope, offset, spline](double input) -> double
            {
                return slope * input + offset + spline->Eval(input);
            };

            calibrations.push_back(calibration_func);
        }

        std::cout << "Loaded " << calibrations.size() << " calibration functions from " << calibrationDir << std::endl;
    }
    catch (const std::exception &e)
    {
        std::cerr << "Error reading calibration directory: " << e.what() << std::endl;
    }

    return calibrations;
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

std::vector<Double_t> readLinearCalParams(const std::string &filename)
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
        double slope, offset;
        if (iss >> slope >> offset)
        {
            params[0] = slope;
            params[1] = offset;
            break;
        }
        // If the line didn't contain two valid numbers, keep searching
    }

    calfile.close();
    return params;
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

    /* #region Calibration Setup */

    // Open calibration file and make calibration functions
    std::vector<std::function<double(double)>> cloverCrossECal(Histograms::kDigitizerChannels), cloverBackECal(Histograms::kDigitizerChannels), posSigECal(Histograms::kDigitizerChannels), cebrAllECal(Histograms::kDigitizerChannels);

    // Load calibration splines and linear params
    std::vector<TSpline3> clover_cross_cal_splines;
    std::vector<std::vector<Double_t>> clover_cross_cal_linear_params;

    for (int det : {1, 3, 5, 7})
    {
        for (int xtal = 1; xtal < 5; xtal++)
        {
            std::string cal_filename = Form("%s/C%iE%i.cal_params.txt", ENERGY_CAL_DIR, det, xtal);
            clover_cross_cal_linear_params.push_back(readLinearCalParams(cal_filename));
            clover_cross_cal_splines.push_back(createCalSpline(cal_filename));
        }
    }

    // Load calibration splines and linear params
    std::vector<TSpline3> clover_back_cal_splines;
    std::vector<std::vector<Double_t>> clover_back_cal_linear_params;
    for (int det : {1, 3, 5, 7})
    {
        for (int xtal = 1; xtal < 5; xtal++)
        {
            std::string cal_filename = Form("%s/B%iE%i.cal_params.txt", ENERGY_CAL_DIR, det, xtal);
            clover_back_cal_linear_params.push_back(readLinearCalParams(cal_filename));
            clover_back_cal_splines.push_back(createCalSpline(cal_filename));
        }
    }

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

        // Clover Back Histograms
        auto cb_amp_raw_hist = clover_back_amp_raw.MakePtr();
        auto cb_cht_raw_hist = clover_back_cht_raw.MakePtr();
        auto cb_plu_hist = clover_back_plu.MakePtr();
        auto cb_trt_raw_hist = clover_back_trt_raw.MakePtr();
        auto cb_mdt_raw_hist = clover_back_mdt_raw.MakePtr();

        auto cb_E_hist = clover_back_E.MakePtr();
        auto cb_cht_hist = clover_back_cht.MakePtr();
        auto cb_mdt_hist = clover_back_mdt.MakePtr();
        auto cb_trt_hist = clover_back_trt.MakePtr();

        auto cb_sum_hist = clover_back_sum.MakePtr();
        auto cb_adb_hist = clover_back_addback.MakePtr();
        auto cb_adb_mult_hist = clover_back_addback_mult.MakePtr();

        // Positive Signal Histograms
        auto ps_amp_raw_hist = pos_sig_amp_raw.MakePtr();
        auto ps_cht_raw_hist = pos_sig_cht_raw.MakePtr();
        auto ps_plu_hist = pos_sig_plu.MakePtr();
        auto ps_trt_raw_hist = pos_sig_trt_raw.MakePtr();
        auto ps_mdt_raw_hist = pos_sig_mdt_raw.MakePtr();

        auto ps_E_hist = pos_sig_E.MakePtr();
        auto ps_cht_hist = pos_sig_cht.MakePtr();
        auto ps_mdt_hist = pos_sig_mdt.MakePtr();
        auto ps_trt_hist = pos_sig_trt.MakePtr();

        auto ps_sum_hist = pos_sig_sum.MakePtr();
        auto ps_adb_hist = pos_sig_addback.MakePtr();
        auto ps_adb_mult_hist = pos_sig_addback_mult.MakePtr();

        // CeBr All Histograms
        auto ce_inl_raw_hist = cebr_all_inl_raw.MakePtr();
        auto ce_cht_raw_hist = cebr_all_cht_raw.MakePtr();
        auto ce_ins_raw_hist = cebr_all_ins_raw.MakePtr();
        auto ce_trt_raw_hist = cebr_all_trt_raw.MakePtr();
        auto ce_mdt_raw_hist = cebr_all_mdt_raw.MakePtr();

        auto ce_El_hist = cebr_all_El.MakePtr();
        auto ce_Es_hist = cebr_all_Es.MakePtr();
        auto ce_cht_hist = cebr_all_cht.MakePtr();
        auto ce_mdt_hist = cebr_all_mdt.MakePtr();
        auto ce_trt_hist = cebr_all_trt.MakePtr();

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
                    double energy = cc_amp[ch]; // cloverCrossECal[ch](cc_amp[ch]);
                    double cht = cc_cht[ch] * kNsPerBin;

                    cc_E_hist->Fill(energy, ch);
                    cc_cht_hist->Fill(cht, ch);
                    cc_sum_hist->Fill(energy, det); // ch / 4 is the detector number
                    xtal_energies.push_back(energy);
                    xtal_times.push_back(cht);
                }

                if (!xtal_energies.empty())
                {
                    cc_adb_hist->Fill(cloverAddBackEnergy(xtal_energies, xtal_times), det);
                    cc_adb_mult_hist->Fill(xtal_energies.size(), det);
                }
            }

            /* #endregion */

            /* #region clover_back */

            // Module Time
            cb_mdt_raw_hist->Fill(cb_mdt[0]);
            cb_mdt_hist->Fill(cb_mdt[0] * kNsPerBin);

            // Trigger Times
            cb_trt_raw_hist->Fill(cb_trt[0], 0);
            cb_trt_hist->Fill(cb_trt[0] * kNsPerBin, 0);
            cb_trt_raw_hist->Fill(cb_trt[1], 1);
            cb_trt_hist->Fill(cb_trt[1] * kNsPerBin, 1);

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
                    cb_amp_raw_hist->Fill(cb_amp[ch], ch);
                    cb_cht_raw_hist->Fill(cb_cht[ch], ch);
                    cb_plu_hist->Fill(cb_plu[ch], ch);

                    // Calibrated Histograms
                    double energy = cb_amp[ch]; // cloverCrossECal[ch](cb_amp[ch]);
                    double cht = cb_cht[ch] * kNsPerBin;

                    cb_E_hist->Fill(energy, ch);
                    cb_cht_hist->Fill(cht, ch);
                    cb_sum_hist->Fill(energy, det); // ch / 4 is the detector number
                    xtal_energies.push_back(energy);
                    xtal_times.push_back(cht);
                }

                if (!xtal_energies.empty())
                {
                    cb_adb_hist->Fill(cloverAddBackEnergy(xtal_energies, xtal_times), det);
                    cb_adb_mult_hist->Fill(xtal_energies.size(), det);
                }
            }

            /* #endregion */

            /* #region pos_sig */

            // Module Time
            ps_mdt_raw_hist->Fill(ps_mdt[0]);
            ps_mdt_hist->Fill(ps_mdt[0] * kNsPerBin);

            // Trigger Times
            ps_trt_raw_hist->Fill(ps_trt[0], 0);
            ps_trt_hist->Fill(ps_trt[0] * kNsPerBin, 0);
            ps_trt_raw_hist->Fill(ps_trt[1], 1);
            ps_trt_hist->Fill(ps_trt[1] * kNsPerBin, 1);

            // Main Loop

            // Detector Loop

            // ZDEG
            ps_amp_raw_hist->Fill(ps_amp[0], 0);
            ps_cht_raw_hist->Fill(ps_cht[0], 0);
            ps_plu_hist->Fill(ps_plu[0], 0);

            // S4E1
            ps_amp_raw_hist->Fill(ps_amp[2], 2);
            ps_cht_raw_hist->Fill(ps_cht[2], 2);
            ps_plu_hist->Fill(ps_plu[2], 2);

            // B4E1-4
            size_t det = 1; // B4 Detector
            std::vector<double> xtal_energies;
            std::vector<double> xtal_times;

            // Crystal Loop

            for (size_t xtal = 0; xtal < 4; xtal++)
            {
                auto ch = det * 4 + xtal; // Channel number 4-7

                // Raw Histograms
                ps_amp_raw_hist->Fill(ps_amp[ch], ch);
                ps_cht_raw_hist->Fill(ps_cht[ch], ch);
                ps_plu_hist->Fill(ps_plu[ch], ch);

                // Calibrated Histograms
                double energy = ps_amp[ch]; // cloverCrossECal[ch](cc_amp[ch]);
                double cht = ps_cht[ch] * kNsPerBin;

                ps_E_hist->Fill(energy, ch);
                ps_cht_hist->Fill(cht, ch);
                ps_sum_hist->Fill(energy, det); // ch / 4 is the detector number
                xtal_energies.push_back(energy);
                xtal_times.push_back(cht);
            }

            if (!xtal_energies.empty())
            {
                ps_adb_hist->Fill(cloverAddBackEnergy(xtal_energies, xtal_times), det);
                ps_adb_mult_hist->Fill(xtal_energies.size(), det);
            }

            /* #endregion */

            /* #region CeBr All Module */

            for (int ch : {cB, cC, cD, cF, cG, cH, cK, cO, cBJ, cBK, cBL, L3, MPAD})
            {
                if (!std::isnan(ce_inl[ch]) && !std::isnan(ce_cht[ch]) && !std::isnan(ce_ins[ch] && !std::isnan(ce_mdt[0])))
                {
                    ce_inl_raw_hist->Fill(ce_inl[ch], ch);
                    ce_cht_raw_hist->Fill(ce_cht[ch], ch);
                    ce_ins_raw_hist->Fill(ce_ins[ch], ch);

                    double energy = ce_inl[ch]; // cebrAllECal[channel](ce_inl[channel]);
                    double ch_time = ce_cht[ch] * kNsPerBin;
                    double md_time = ce_mdt[0] * kNsPerBin;

                    if (energy > CEBR_ALL_ENERGY_THRESHOLD)
                    {
                        ce_El_hist->Fill(energy, ch);
                        ce_cht_hist->Fill(ch_time, ch);
                        ce_mdt_hist->Fill(md_time);
                    }
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
    clover_cross_cht_raw.Write();
    clover_cross_plu.Write();
    clover_cross_trt_raw.Write();
    clover_cross_mdt_raw.Write();

    clover_cross_E.Write();
    clover_cross_cht.Write();
    clover_cross_mdt.Write();
    clover_cross_trt.Write();

    clover_cross_sum.Write();
    clover_cross_addback.Write();
    clover_cross_addback_mult.Write();

    outfile->cd();

    // Clover Back Histograms

    auto cb_dir = outfile->mkdir("clover_back");
    cb_dir->cd();
    clover_back_amp_raw.Write();
    clover_back_cht_raw.Write();
    clover_back_plu.Write();
    clover_back_trt_raw.Write();
    clover_back_mdt_raw.Write();

    clover_back_E.Write();
    clover_back_cht.Write();
    clover_back_mdt.Write();
    clover_back_trt.Write();
    clover_back_sum.Write();
    clover_back_addback.Write();
    clover_back_addback_mult.Write();

    outfile->cd();

    // Positive Signal Histograms
    auto ps_dir = outfile->mkdir("pos_sig");
    ps_dir->cd();
    pos_sig_amp_raw.Write();
    pos_sig_cht_raw.Write();
    pos_sig_plu.Write();
    pos_sig_trt_raw.Write();
    pos_sig_mdt_raw.Write();

    pos_sig_E.Write();
    pos_sig_cht.Write();
    pos_sig_mdt.Write();
    pos_sig_trt.Write();
    pos_sig_sum.Write();
    pos_sig_addback.Write();
    pos_sig_addback_mult.Write();

    outfile->cd();

    // CeBr All Histograms
    auto ce_dir = outfile->mkdir("cebr_all");
    ce_dir->cd();
    cebr_all_inl_raw.Write();
    cebr_all_cht_raw.Write();
    cebr_all_ins_raw.Write();
    cebr_all_trt_raw.Write();
    cebr_all_mdt_raw.Write();

    cebr_all_El.Write();
    cebr_all_Es.Write();
    cebr_all_cht.Write();
    cebr_all_mdt.Write();
    cebr_all_trt.Write();

    outfile->cd();

    /* #endregion */

    std::cout << "Saved histograms to file: " << outfile->GetName() << std::endl;

    outfile->Close();
    delete outfile;

    std::cout << "Done!" << std::endl;

    return 0;
}