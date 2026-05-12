// C++ Includes
#include <chrono>
#include <fstream>
#include <iostream>
#include <sstream>
#include <thread>
#include <vector>

// Boost Includes
#include <boost/program_options.hpp>

// ROOT Includes
#include <TString.h>

// Project Includes
#include "CAConfiguration.hpp"
#include "CAUtilities.hpp"

namespace po = boost::program_options;

CAUtilities::Args CAUtilities::ParseArguments(int argc, char* argv[])
{
    po::options_description desc("Options");
    desc.add_options()
        ("help,h", "Display this help message")
        ("mode,m", po::value<std::string>()->default_value("raw"), "Processing mode (raw, cal, xtcorr)")
        ("caldir,c", po::value<std::string>(), "Directory containing calibration files [cal, xtcorr only]")
        ("gsfile,g", po::value<std::string>(), "File containing gain shift data [cal, xtcorr only]")
        ("xtfile,x", po::value<std::string>(), "File containing crosstalk correction matrices [xtcorr only]")
        ("run-file,r", po::value<std::string>()->required(), "Input run file name")
        ("out-file,o", po::value<std::string>()->required(), "Output file name");

    po::positional_options_description pos;
    pos.add("run-file", 1);
    pos.add("out-file", 1);

    po::variables_map vm;
    try
    {
        po::store(po::command_line_parser(argc, argv).options(desc).positional(pos).run(), vm);

        if (vm.count("help"))
        {
            printf("Usage: %s [options] -r <run_file_name> -o <output_file_name>\n\n", argv[0]);
            std::cout << desc << std::endl;
            exit(EXIT_SUCCESS);
        }

        po::notify(vm);

        // Validate mode-dependent options
        const std::string mode = vm["mode"].as<std::string>();
        const bool needsCalibration = (mode == "cal" || mode == "xtcorr");
        if (needsCalibration && !vm.count("gsfile"))
            throw po::error("-g [ --gsfile ] is required for mode '" + mode + "'");
        if (needsCalibration && !vm.count("caldir"))
            throw po::error("-c [ --caldir ] is required for mode '" + mode + "'");
        if (!needsCalibration && (vm.count("gsfile") || !vm["caldir"].defaulted()))
            std::cerr << "[WARN] --caldir and --gsfile are ignored in mode '" << mode << "'\n";
    }
    catch (const po::error& e)
    {
        std::cerr << "[ERROR] " << e.what() << "\n\n";
        printf("Usage: %s [options] -r <run_file_name> -o <output_file_name>\n\n", argv[0]);
        std::cout << desc << std::endl;
        exit(EXIT_FAILURE);
    }

    Args args;
    args.mode = vm["mode"].as<std::string>();
    args.calibrationDir = vm.count("caldir") ? vm["caldir"].as<std::string>() : "";
    args.gainShiftFile = vm.count("gsfile") ? vm["gsfile"].as<std::string>() : "";
    args.xtalkFile = vm.count("xtfile") ? vm["xtfile"].as<std::string>() : "";
    args.runFileName = vm["run-file"].as<std::string>();
    args.outputFileName = vm["out-file"].as<std::string>();

    return args;
}

void CAUtilities::PrintConfiguration(const Args& args)
{
    const bool needsCalibration = (args.mode == "cal" || args.mode == "xtcorr");

    std::cout << "--------------- Current Configuration ------------------" << std::endl;
    std::cout << "Processing mode: " << args.mode << std::endl;
    std::cout << "Run file:        " << args.runFileName << std::endl;
    std::cout << "Output file:     " << args.outputFileName << std::endl;
    if (needsCalibration)
    {
        std::cout << "Calibration dir: " << args.calibrationDir << std::endl;
        std::cout << "Gain-shift file: " << args.gainShiftFile << std::endl;
    }
    std::cout << "Max Threads:     " << kMaxThreads << std::endl;
    std::cout << "--------------------------------------------------------" << std::endl;
}

void CAUtilities::DisplayProgressBar(std::atomic<uint64_t>& processedEntries, uint64_t totalEntries)
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

std::vector<std::vector<std::vector<double>>> CAUtilities::ReadCAFile(const std::string& fileName)
{
    std::vector<std::vector<std::vector<double>>> data;
    std::ifstream inputFile(fileName);
    if (!inputFile.is_open())
    {
        throw std::runtime_error("[ERROR] Could not open file " + fileName);
    }

    std::string line;
    std::string currentSection;
    while (std::getline(inputFile, line))
    {
        // Skip empty lines
        if (line.empty())
        {
            continue;
        }

        int currentSectionIndex = -1;
        // Check for comment/section header lines
        if (line[0] == '#')
        {
            // Check if it's a section header (not the column header)
            if (line.find("Channel") == std::string::npos)
            {
                currentSection = line.substr(2); // Remove "# "
                #if DEBUG >= 2
                printf("Reading section: %s\n", currentSection.c_str());
                #endif
                data.push_back(std::vector<std::vector<double>>());
            }
            continue;
        }

        // Parse data line: Channel, Val1, Val2, ...
        data.back().push_back(std::vector<double>());
        std::istringstream iss(line);
        int channel;
        double value;
        while (iss >> channel)
        {
            data.back().back().push_back(channel);
            while (iss >> value)
            {
                data.back().back().push_back(value);
            }
        }
    }
    inputFile.close();
    return data;
}