#ifndef CAADDBACK_HPP
#define CAADDBACK_HPP

// C++ Includes
#include <array>

// ROOT Includes

// Project Includes

namespace CAAddBack
{
    // Constants
    static constexpr double kAddBackThreshold = 0; // Energy threshold (keV) for add-back
    static constexpr double kAddBackWindow = 150;  // Time window (ns) around primary hit for add-back

    double GetAddBackEnergy(std::array<double, 4> xtalE, std::array<double, 4> xtalT);

} // namespace CAAddBack

#endif // CAADDBACK_HPP