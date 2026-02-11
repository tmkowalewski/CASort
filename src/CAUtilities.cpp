// C++ Includes
#include <chrono>
#include <fstream>
#include <iostream>
#include <sstream>
#include <thread>
#include <vector>

// ROOT Includes
#include <TString.h>

// Project Includes
#include "CAConfiguration.hpp"
#include "CAUtilities.hpp"

CAUtilities::Args CAUtilities::ParseArguments(int argc, char* argv[])
{
    if (argc != 6)
    {
        printf("Usage: %s <calibration_dir> <gain_shift_dir> <run_file_dir> <run_number> <output_file_name>\n", argv[0]);
        exit(EXIT_FAILURE);
    }

    Args args;
    args.calibrationDir = argv[1];
    args.gainShiftDir = argv[2];
    std::string runFileDir = argv[3];
    args.runNumber = std::stoi(argv[4]);
    std::string partialRunFileName = Form(RUN_FILE_NAME_TEMPLATE, args.runNumber);
    args.runFileName = Form("%s/%s", runFileDir.c_str(), partialRunFileName.c_str());
    args.outputFileName = argv[5];

    return args;
}

void CAUtilities::PrintConfiguration(const Args& args)
{
    std::cout << "--------------- Current Configuration ------------------" << std::endl;
    std::cout << "Calibration directory: " << args.calibrationDir << std::endl;
    std::cout << "Gain-shift directory: " << args.gainShiftDir << std::endl;
    std::cout << "Run file: " << args.runFileName << std::endl;
    std::cout << "Output file: " << args.outputFileName << std::endl;
    std::cout << "Max Threads: " << kMaxThreads << std::endl;
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
        std::cerr << "Error: Could not open file " << fileName << std::endl;
        return data;
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