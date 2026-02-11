#ifndef CAUTILITIES_HPP
#define CAUTILITIES_HPP

// C++ Includes
#include <atomic>
#include <string>
#include <vector>

// ROOT Includes

// Project Includes

namespace CAUtilities
{
    struct Args
    {
        std::string calibrationDir;
        std::string gainShiftDir;
        std::string runFileName;
        std::string outputFileName;
        int runNumber;
    };

    Args ParseArguments(int argc, char* argv[]);

    void PrintConfiguration(const Args& args);

    void DisplayProgressBar(std::atomic<uint64_t>& processedEntries, uint64_t totalEntries);

    std::vector<std::vector<std::vector<double>>> ReadCAFile(const std::string& fileName);

} // namespace CAUtilities

#endif // CAUTILITIES_HPP