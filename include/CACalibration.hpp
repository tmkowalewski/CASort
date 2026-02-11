#ifndef CACALIBRATION_HPP
#define CACALIBRATION_HPP

// C++ Includes
#include <functional>
#include <string>
#include <vector>

// ROOT Includes
#include <TSpline.h>

// Project Includes

namespace CACalibration
{
    inline constexpr double kMaxCalibrationEnergy = 7282.92; // Maximum energy for calibration spline

    TSpline3 LoadSplineCorrParams(const std::string& fileName);

    std::vector<double> LoadLinearCalParams(const std::string& fileName);

    std::function<double(double)> MakeCalibration(const std::string& fileName);

} // namespace CACalibration

#endif // CACALIBRATION_HPP