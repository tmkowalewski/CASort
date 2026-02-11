#ifndef CAGAINCORRECTION_HPP
#define CAGAINCORRECTION_HPP

// C++ Includes
#include <functional>
#include <string>
#include <vector>

// ROOT Includes

// Project Includes

namespace CAGainCorrection
{

    inline std::string gGainCorrectionDir = "";

    std::vector<std::vector<std::function<double(double)>>> MakeCorrections(const std::string& gainshiftDir, const unsigned int runNumber);

} // namespace CAGainCorrection

#endif // CAGAINCORRECTION_HPP