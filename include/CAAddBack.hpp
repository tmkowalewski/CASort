#ifndef CAADDBACK_HPP
#define CAADDBACK_HPP

// C++ Includes
#include <array>

// ROOT Includes

// Project Includes

namespace CAAddBack
{
    // Constants
    static constexpr double kAddBackThreshold = 150; // Energy threshold (keV) for add-back
    static constexpr double kAddBackWindow = 150;    // Time window (ns) around primary hit for add-back

    double GetAddBackEnergy(std::array<double, 4> xtal_E, std::array<double, 4> xtal_T);

} // namespace CAAddBack

#endif // CAADDBACK_HPP