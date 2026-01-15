// C++ Includes
#include <iostream>
#include <thread>
#include <chrono>

// ROOT Includes

// Project Includes
#include "Utilities.hpp"

void Utilities::DisplayProgressBar(std::atomic<uint64_t> &processedEntries, uint64_t totalEntries)
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