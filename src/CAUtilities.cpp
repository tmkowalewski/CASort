// C++ Includes
#include <iostream>
#include <fstream>
#include <sstream>
#include <thread>
#include <chrono>
#include <vector>

// ROOT Includes

// Project Includes
#include "CAConfiguration.hpp"
#include "CAUtilities.hpp"

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

using vector3D = std::vector<std::vector<std::vector<double>>>;
vector3D CAUtilities::ReadCAFile(const std::string& filename)
{
    vector3D data;
    std::ifstream infile(filename);
    if (!infile.is_open())
    {
        std::cerr << "Error: Could not open file " << filename << std::endl;
        return data;
    }

    std::string line;
    std::string current_section;
    while (std::getline(infile, line))
    {
        // Skip empty lines
        if (line.empty())
        {
            continue;
        }

        int current_section_index = -1;
        // Check for comment/section header lines
        if (line[0] == '#')
        {
            // Check if it's a section header (not the column header)
            if (line.find("Channel") == std::string::npos)
            {
                current_section = line.substr(2); // Remove "# "
                #if DEBUG >= 2
                printf("Reading section: %s\n", current_section.c_str());
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
        while (iss >> channel >> value)
        {
            data.back().back().push_back(channel);
            data.back().back().push_back(value);
            while (iss >> value)
            {
                data.back().back().push_back(value);
            }
        }
    }
    infile.close();
    return data;
}