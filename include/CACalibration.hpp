#ifndef CACALIBRATION_HPP
#define CACALIBRATION_HPP

// C++ Includes
#include <string>
#include <vector>
#include <functional>

// ROOT Includes
#include <TSpline.h>

// Project Includes

namespace CACalibration
{
    inline constexpr double kMaxCalibrationEnergy = 7282.92; // Maximum energy for calibration spline

    TSpline3 CreateSplineCorrection(const std::string& filename);

    std::vector<double> ReadLinearCalParams(const std::string& filename);

    std::function<double(double)> MakeCalibration(const std::string& filename);

} // namespace CACalibration

#endif // CACALIBRATION_HPP