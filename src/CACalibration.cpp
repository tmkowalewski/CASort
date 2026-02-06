// C++ Includes
#include <fstream>
#include <iostream>
#include <sstream>

// ROOT Includes

// Project Includes
#include "CACalibration.hpp"

TSpline3 CACalibration::LoadSplineCorrParams(const std::string& filename)
{
    std::vector<double> knot_x, knot_y;
    std::ifstream calfile(filename);
    if (!calfile.is_open())
    {
        std::cerr << "Failed to open calibration file: " << filename << std::endl;
        return TSpline3();
    }

    // Read spline knots from the calibration file
    std::string line;
    while (std::getline(calfile, line))
    {
        bool skipped_linear = false;
        // Helper lambda to trim whitespace
        auto trim = [](std::string& s)
            {
                const char* whitespace = " \t\n\r\f\v";
                size_t start = s.find_first_not_of(whitespace);
                if (start == std::string::npos)
                {
                    s.clear();
                    return;
                }
                size_t end = s.find_last_not_of(whitespace);
                s = s.substr(start, end - start + 1);
            };

        while (std::getline(calfile, line))
        {
            trim(line);
            if (line.empty())
                continue;
            if (line[0] == '#')
                continue;

            // The first non-comment/non-empty line contains the linear parameters (slope offset).
            // Skip that line and only start collecting knots afterwards.
            if (!skipped_linear)
            {
                skipped_linear = true;
                continue;
            }

            std::istringstream iss(line);
            double x, y;
            if (iss >> x >> y)
            {
                knot_x.push_back(x);
                knot_y.push_back(y);
            }
        }
        if (line.empty() || line[0] == '#')
            continue;

        std::istringstream iss(line);
        double x, y;
        if (iss >> x >> y)
        {
            knot_x.push_back(x);
            knot_y.push_back(y);
        }
    }

    calfile.close();

    if (knot_x.empty())
    {
        std::cerr << "No spline knots found in " << filename << std::endl;
        return TSpline3();
    }

    // Create and return the spline
    return TSpline3("spline", knot_x.data(), knot_y.data(), knot_x.size(), "b1e1");
}

std::vector<double> CACalibration::LoadLinearCalParams(const std::string& filename)
{
    std::vector<Double_t> params(2, 0.0); // Initialize with two zeros
    std::ifstream calfile(filename);
    if (!calfile.is_open())
    {
        std::cerr << "Failed to open calibration file: " << filename << std::endl;
        return params;
    }

    std::string line;
    while (std::getline(calfile, line))
    {
        // Trim leading whitespace
        size_t first = line.find_first_not_of(" \t\r\n");
        if (first == std::string::npos)
            continue; // empty line

        // Skip comment lines
        if (line[first] == '#')
            continue;

        std::istringstream iss(line);
        double offset, slope;
        if (iss >> offset >> slope)
        {
            params[0] = offset;
            params[1] = slope;
            break;
        }
        // If the line didn't contain two valid numbers, keep searching
    }

    calfile.close();
    return params;
}

std::function<double(double)> CACalibration::MakeCalibration(const std::string& filename)
{
    auto linear_params = LoadLinearCalParams(filename);
    auto cal_spline = LoadSplineCorrParams(filename);

    double offset = linear_params[0];
    double slope = linear_params[1];

    // Create calibration function: output = slope * input + offset + spline(input)
    auto calibration_func = [slope, offset, cal_spline](double input) -> double
        {
            double lincal_E = slope * input + offset;
            double spline_corr = cal_spline.Eval(lincal_E);
            double energy = input < kMaxCalibrationEnergy ? lincal_E + spline_corr : lincal_E; // Only trust the spline when interpolating
            return energy;
        };

    // Return the callable calibration lambda instead of an empty std::function
    return calibration_func;
}