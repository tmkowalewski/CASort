#ifndef CAUTILITIES_HPP
#define CAUTILITIES_HPP

// C++ Includes
#include <atomic>
#include <vector>
#include <string>

// ROOT Includes

// Project Includes

namespace CAUtilities
{
    void DisplayProgressBar(std::atomic<uint64_t>& processedEntries, uint64_t totalEntries);

    std::vector<std::vector<std::vector<double>>> ReadCAFile(const std::string& filename);

} // namespace CAUtilities

#endif // CAUTILITIES_HPP