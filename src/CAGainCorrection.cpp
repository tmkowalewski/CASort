// C++ Includes
#include <filesystem>
#include <fstream>
#include <iostream>
#include <sstream>

// ROOT Includes
#include <TString.h>

// Project Includes
#include "CAGainCorrection.hpp"
#include "CAUtilities.hpp"

std::vector<std::vector<std::function<double(double)>>> CAGainCorrection::MakeCorrections(const std::string& fileName)
{
    std::vector<std::vector<std::function<double(double)>>> gainshiftFunctions;

    std::string gsFileName = fileName;

    if (fileName.empty())
    {
        printf("[WARN] Gain shift file \"%s\" not found! Using default gain shift values.\n", fileName.c_str());
        gsFileName = Form("70Ge_default.cags");
    }
    printf("[INFO] Loading gain shift data from file %s\n", gsFileName.c_str());
    auto rawData = CAUtilities::ReadCAFile(gsFileName);
    for (const auto& module_Data : rawData)
    {
        std::vector<std::function<double(double)>> channelFunctions;
        for (const auto& channelData : module_Data)
        {

            if (channelData.size() != 3)
            {
                std::cout << channelData.size() << std::endl;
                throw std::runtime_error(Form("[ERROR] Unexpected data format in gain shift file %s. Correct format:\n channel_number offset gain\n", fileName.c_str()));
            }
            int channel = static_cast<int>(channelData[0]);
            double offset = channelData[1];
            double gain = channelData[2];
            channelFunctions.push_back([offset, gain](double x)
                                       { return gain * x + offset; });
#if DEBUG >= 2
            printf("[INFO] Channel %d: Offset = %.6f, Gain = %.6f\n", channel, offset, gain);
#endif
        }
        gainshiftFunctions.push_back(channelFunctions);
    }

    return gainshiftFunctions;
}