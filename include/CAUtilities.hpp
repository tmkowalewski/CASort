#ifndef CAUTILITIES_HPP
#define CAUTILITIES_HPP

// C++ Includes
#include <atomic>

// ROOT Includes

// Project Includes

namespace CAUtilities
{
    void DisplayProgressBar(std::atomic<uint64_t> &processedEntries, uint64_t totalEntries);
} // namespace CAUtilities

#endif // CAUTILITIES_HPP