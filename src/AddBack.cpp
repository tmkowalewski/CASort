// C++ Includes
#include <iostream>
#include <cmath>

// ROOT Includes

// Project Includes
#include "AddBack.hpp"

double AddBack::GetAddBackEnergy(std::vector<double> crystal_energies, std::vector<double> crystal_times)
{
    // Clover add-back function, modified from Samantha's code

    if (crystal_energies.size() != crystal_times.size())
    {
        std::cerr << "Error: GetAddBackEnergy requires the same number of detector energies and times" << std::endl;
        return 0;
    }

    double highest_energy = 0;
    int highest_energy_index = -1;
    double final_energy = 0;
    std::vector<double> delta_t;
    unsigned short mult = 0; // Addback multiplicity counter

    // Find highest energy hit, this is most likely the primary hit
    for (size_t crystal = 0; crystal < crystal_energies.size(); crystal++)
    {
        if (crystal_energies[crystal] > kAddBackThreshold) // Energy cut in keV
        {
            if (crystal_energies[crystal] > highest_energy)
            {
                highest_energy = crystal_energies[crystal];
                highest_energy_index = crystal;
            }
        }
    }
    if (highest_energy_index > -1) // If a primary hit was found
    {
        for (size_t crystal = 0; crystal < crystal_energies.size(); crystal++) // Perform the addback
        {
            delta_t.push_back(fabs((crystal_times[highest_energy_index] - crystal_times[crystal])));
            if (fabs((crystal_times[highest_energy_index] - crystal_times[crystal])) < kAddBackWindow) // Time cut in ns
            {
                if (crystal_energies[crystal] > kAddBackThreshold) // Energy cut in keV
                {
                    final_energy += crystal_energies[crystal];
                    mult++;
                }
            }
        }
    }

    // fmt::print("Energies: {}\n Times: {}\n Final Energy: {}\n\n", crystal_energies, delta_t, final_energy);

    return final_energy;
}