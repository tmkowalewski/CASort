#ifndef CALIBRATION_HPP
#define CALIBRATION_HPP

// C++ Includes
#include <string>
#include <vector>
#include <functional>

// ROOT Includes
#include <TSpline.h>

// Project Includes

namespace Calibration
{
    TSpline3 CreateSplineCorrection(const std::string &filename);

    std::vector<double> ReadLinearCalParams(const std::string &filename);

    std::function<double(double)> MakeCalibration(std::vector<double> linear_params, TSpline3 cal_spline);

} // namespace Calibration

#endif // CALIBRATION_HPP