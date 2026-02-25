// C++ Includes
#include <cmath>
#include <iostream>

// ROOT Includes

// Project Includes
#include "CAAddBack.hpp"

double CAAddBack::GetAddBackEnergy(std::array<double, 4> xtalE, std::array<double, 4> xtalT)
{
    // Clover add-back function, modified from Samantha's code

    double primaryE = kAddBackThreshold; // Start with threshold, if no hits are above this then add-back will return 0
    int primaryIdx = -1;
    double finalE = 0;
    std::array<double, 4> deltaT;

    // Find highest energy hit, this is most likely the primary hit
    for (size_t xtal = 0; xtal < 4; xtal++)
    {
        if (xtalE[xtal] >= primaryE)
        {
            primaryE = xtalE[xtal];
            primaryIdx = xtal;
        }
    }
    if (primaryIdx > -1) // If a primary hit was found
    {
        const double primaryTime = xtalT[primaryIdx];
        finalE = primaryE;                      // Primary always passes both checks, add it directly
        for (size_t xtal = 0; xtal < 4; xtal++) // Perform the addback
        {
            if (xtal == static_cast<size_t>(primaryIdx))
                continue;
            if (xtalE[xtal] > kAddBackThreshold &&
                fabs(primaryTime - xtalT[xtal]) < kAddBackWindow)
            {
                finalE += xtalE[xtal];
            }
        }
    }
#if DEBUG > 1
    else
    {

        std::cout << "[WARNING] No primary hit found for add-back. Returning 0 energy." << std::endl;
        std::cout << "xtalE = {";
        for (size_t i = 0; i < xtalE.size(); i++)
            std::cout << xtalE[i] << ", ";
        std::cout << "}" << std::endl;
    }
#endif // DEBUG

    return finalE;
}