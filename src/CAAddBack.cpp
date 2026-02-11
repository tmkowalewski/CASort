// C++ Includes
#include <cmath>
#include <iostream>

// ROOT Includes

// Project Includes
#include "CAAddBack.hpp"

double CAAddBack::GetAddBackEnergy(std::array<double, 4> xtalE, std::array<double, 4> xtalT)
{
    // Clover add-back function, modified from Samantha's code

    double primaryE = 0;
    int primaryIdx = -1;
    double finalE = 0;
    std::array<double, 4> deltaT;

    // Find highest energy hit, this is most likely the primary hit
    for (size_t xtal = 0; xtal < xtalE.size(); xtal++)
    {
        if (xtalE[xtal] > kAddBackThreshold) // Energy cut in keV
        {
            if (xtalE[xtal] > primaryE)
            {
                primaryE = xtalE[xtal];
                primaryIdx = xtal;
            }
        }
    }
    if (primaryIdx > -1) // If a primary hit was found
    {
        for (size_t xtal = 0; xtal < xtalE.size(); xtal++) // Perform the addback
        {
            deltaT[xtal] = fabs((xtalT[primaryIdx] - xtalT[xtal]));
            if (fabs((xtalT[primaryIdx] - xtalT[xtal])) < kAddBackWindow && xtalE[xtal] > kAddBackThreshold) // Time cut in ns, energy cut in keV
            {
                if (xtalE[xtal] > kAddBackThreshold) // Energy cut in keV
                {
                    finalE += xtalE[xtal];
                }
            }
        }
    }

    return finalE;
}