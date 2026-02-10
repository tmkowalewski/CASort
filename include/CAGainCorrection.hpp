#ifndef CAGAINCORRECTION_HPP
#define CAGAINCORRECTION_HPP

// C++ Includes
#include <string>
#include <vector>
#include <functional>

// ROOT Includes

// Project Includes

namespace CAGainCorrection
{

    inline std::string gGainCorrectionDir = "";

    std::vector<std::vector<std::function<double(double)>>> MakeCorrections(const std::string &gain_shift_dir, const unsigned int run_number);

} // namespace CAGainCorrection

#endif // CAGAINCORRECTION_HPP