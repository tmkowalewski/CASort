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
        const char* calibrationDir;
        const char* gainShiftDir;
        const char* runFileName;
        const char* outputFileName;
        int runNumber;
    };

    Args ParseArguments(int argc, char* argv[]);

    void PrintConfiguration(const Args& args);

    void DisplayProgressBar(std::atomic<uint64_t>& processedEntries, uint64_t totalEntries);

    std::vector<std::vector<std::vector<double>>> ReadCAFile(const std::string& filename);

} // namespace CAUtilities

#endif // CAUTILITIES_HPP