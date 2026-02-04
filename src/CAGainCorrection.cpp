// C++ Includes
#include <iostream>
#include <fstream>
#include <sstream>
#include <filesystem>

// ROOT Includes
#include <TString.h>


// Project Includes
#include "CAGainCorrection.hpp"
#include "CAUtilities.hpp"

std::vector<std::vector<std::function<double(double)>>> CAGainCorrection::MakeCorrections(const std::string& gain_shift_dir, const unsigned int run_number)
{
    std::vector<std::vector<std::function<double(double)>>> gain_shift_data;

    // Construct the file path
    std::string run_pattern = Form("run%03d", run_number);
    std::string filename;

    // Find file containing the run pattern in the directory
    for (const auto& entry : std::filesystem::directory_iterator(gain_shift_dir))
    {
        if (entry.path().filename().string().find(run_pattern) != std::string::npos)
        {
            filename = entry.path().string();
            break;
        }
    }
    if (filename.empty())
    {
        throw std::runtime_error(Form("[ERROR] Gain shift file for run %03d not found in directory %s", run_number, gain_shift_dir.c_str()));
    }
    printf("[INFO] Loading gain shift data from file %s\n", filename.c_str());
    auto raw_data = CAUtilities::ReadCAFile(filename);
    for (const auto& module_Data : raw_data)
    {
        std::vector<std::function<double(double)>> channel_funcs;
        for (const auto& channel_data : module_Data)
        {
            if (channel_data.size() != 3)
            {
                throw std::runtime_error(Form("[ERROR] Unexpected data format in gain shift file %s. Correct format:\n channel_number offset gain\n", filename.c_str()));
            }
            int channel = static_cast<int>(channel_data[0]);
            double offset = channel_data[1];
            double gain = channel_data[2];
            channel_funcs.push_back([offset, gain](double x)
                {
                    return gain * x + offset;
                });
            #if DEBUG >= 2
            printf("[INFO] Channel %d: Offset = %.6f, Gain = %.6f\n", channel, offset, gain);
            #endif
        }
        gain_shift_data.push_back(channel_funcs);
    }


    return gain_shift_data;
}