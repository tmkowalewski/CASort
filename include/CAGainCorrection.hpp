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
    std::vector<std::vector<std::function<double(double)>>> GetGainCorrection(const std::string& gain_shift_dir, const unsigned int run_number);

} // namespace CAGainCorrection

#endif // CAGAINCORRECTION_HPP