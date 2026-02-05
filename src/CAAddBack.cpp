// C++ Includes
#include <iostream>
#include <cmath>

// ROOT Includes

// Project Includes
#include "CAAddBack.hpp"

double CAAddBack::GetAddBackEnergy(std::array<double, 4> crystal_energies, std::array<double, 4> xtal_T)
{
    // Clover add-back function, modified from Samantha's code

    double primary_E = 0;
    int primary_idx = -1;
    double final_energy = 0;
    std::array<double, 4> delta_t;
    unsigned short mult = 0; // Addback multiplicity counter

    // Find highest energy hit, this is most likely the primary hit
    for (size_t xtal = 0; xtal < crystal_energies.size(); xtal++)
    {
        if (crystal_energies[xtal] > kAddBackThreshold) // Energy cut in keV
        {
            if (crystal_energies[xtal] > primary_E)
            {
                primary_E = crystal_energies[xtal];
                primary_idx = xtal;
            }
        }
    }
    if (primary_idx > -1) // If a primary hit was found
    {
        for (size_t xtal = 0; xtal < crystal_energies.size(); xtal++) // Perform the addback
        {
            delta_t[xtal] = fabs((xtal_T[primary_idx] - xtal_T[xtal]));
            if (fabs((xtal_T[primary_idx] - xtal_T[xtal])) < kAddBackWindow) // Time cut in ns
            {
                if (crystal_energies[xtal] > kAddBackThreshold) // Energy cut in keV
                {
                    final_energy += crystal_energies[xtal];
                    mult++;
                }
            }
        }
    }

    // fmt::print("Energies: {}\n Times: {}\n Final Energy: {}\n\n", crystal_energies, delta_t, final_energy);

    return final_energy;
}