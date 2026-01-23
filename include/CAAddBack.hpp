#ifndef CAADDBACK_HPP
#define CAADDBACK_HPP

// C++ Includes
#include <vector>

// ROOT Includes

// Project Includes

namespace CAAddBack
{
    // Constants
    static constexpr double kAddBackThreshold = 150; // Energy threshold (keV) for add-back
    static constexpr double kAddBackWindow = 150;    // Time window (ns) around primary hit for add-back

    double GetAddBackEnergy(std::vector<double> crystal_energies, std::vector<double> crystal_times);

} // namespace CAAddBack

#endif // CAADDBACK_HPP