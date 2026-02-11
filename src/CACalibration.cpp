// C++ Includes
#include <fstream>
#include <iostream>
#include <sstream>

// ROOT Includes

// Project Includes
#include "CACalibration.hpp"

TSpline3 CACalibration::LoadSplineCorrParams(const std::string& fileName)
{
    std::vector<double> knotX, knotY;
    std::ifstream calFile(fileName);
    if (!calFile.is_open())
    {
        std::cerr << "Failed to open calibration file: " << fileName << std::endl;
        return TSpline3();
    }

    // Read spline knots from the calibration file
    std::string line;
    while (std::getline(calFile, line))
    {
        bool skippedLinear = false;
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

        while (std::getline(calFile, line))
        {
            trim(line);
            if (line.empty())
                continue;
            if (line[0] == '#')
                continue;

            // The first non-comment/non-empty line contains the linear parameters (slope offset).
            // Skip that line and only start collecting knots afterwards.
            if (!skippedLinear)
            {
                skippedLinear = true;
                continue;
            }

            std::istringstream iss(line);
            double x, y;
            if (iss >> x >> y)
            {
                knotX.push_back(x);
                knotY.push_back(y);
            }
        }
        if (line.empty() || line[0] == '#')
            continue;

        std::istringstream iss(line);
        double x, y;
        if (iss >> x >> y)
        {
            knotX.push_back(x);
            knotY.push_back(y);
        }
    }

    calFile.close();

    if (knotX.empty())
    {
        std::cerr << "No spline knots found in " << fileName << std::endl;
        return TSpline3();
    }

    // Create and return the spline
    return TSpline3("spline", knotX.data(), knotY.data(), knotX.size(), "b1e1");
}

std::vector<double> CACalibration::LoadLinearCalParams(const std::string& fileName)
{
    std::vector<Double_t> params(2, 0.0); // Initialize with two zeros
    std::ifstream calFile(fileName);
    if (!calFile.is_open())
    {
        std::cerr << "Failed to open calibration file: " << fileName << std::endl;
        return params;
    }

    std::string line;
    while (std::getline(calFile, line))
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

    calFile.close();
    return params;
}

std::function<double(double)> CACalibration::MakeCalibration(const std::string& fileName)
{
    auto linearParams = LoadLinearCalParams(fileName);
    auto calSpline = LoadSplineCorrParams(fileName);

    double offset = linearParams[0];
    double slope = linearParams[1];

    // Create calibration function: output = slope * input + offset + spline(input)
    auto calFunc = [slope, offset, calSpline](double input) -> double
    {
        double linearCalE = slope * input + offset;
        double splineCorr = calSpline.Eval(linearCalE);
        double energy = input < kMaxCalibrationEnergy ? linearCalE + splineCorr : linearCalE; // Only trust the spline when interpolating
        return energy;
    };

    // Return the callable calibration lambda instead of an empty std::function
    return calFunc;
}