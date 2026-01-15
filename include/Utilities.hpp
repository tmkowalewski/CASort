#ifndef UTILITIES_HPP
#define UTILITIES_HPP

// C++ Includes
#include <atomic>

// ROOT Includes

// Project Includes

namespace Utilities
{
    void DisplayProgressBar(std::atomic<uint64_t> &processedEntries, uint64_t totalEntries);
} // namespace Utilities

#endif // UTILITIES_HPP