#ifndef CACROSSTALKCORRECTION_HPP
#define CACROSSTALKCORRECTION_HPP

// C++ Includes
#include <vector>

// ROOT Includes
#include <TGraphErrors.h>
#include <TH2D.h>
#include <TMatrixD.h>

// Project Includes
#include "TCAHistogram.hpp"

namespace CACrosstalkCorrection
{
    enum CrystalPair
    {
        E1_E2 = 0,
        E1_E3 = 1,
        E1_E4 = 2,
        E2_E3 = 3,
        E2_E4 = 4,
        E3_E4 = 5
    };

    struct CrosstalkFit
    {
        bool valid = false;
        double alphaXY = std::numeric_limits<double>::quiet_NaN(); // Crosstalk coefficient from x to y
        double alphaYX = std::numeric_limits<double>::quiet_NaN(); // Crosstalk coefficient from y to x
        double alphaXYErr = std::numeric_limits<double>::quiet_NaN();
        double alphaYXErr = std::numeric_limits<double>::quiet_NaN();
        double chi2 = std::numeric_limits<double>::quiet_NaN();
        double ndf = std::numeric_limits<double>::quiet_NaN();
    };

    inline std::string gCrosstalkCorrectionDir = "";

    const double kTargetEnergy = 5018.98;    // Energy in keV of gamma ray used for crosstalk calibration
    const double kEnergyWindow = 10.0;       // Half-width of energy window in keV around target energy
    const double kEnergyCut = 244.0;         // Minimum energy in keV for cut in crosstalk plots
    const double kFitWindow = 30.0;          // Width of fit window in keV around target energy
    const unsigned int kMinCountsPerBin = 1; // Minimum counts per x-bin slice for fitting

    // Function used to model crosstalk effect
    double CrosstalkFitFunction(double* x, double* par);

    void FillXTalkHistograms(const std::array<std::shared_ptr<TH2D>, 6>& xtalkPairHists, const std::array<double, 4>& xtalE, std::array<double, 4>& xtalT);

    std::shared_ptr<TGraphErrors> BuildCrosstalkGraph(const TH2D* hist);

    CrosstalkFit FitCrosstalkCorrection(const TH2D* hist);

    TMatrixD BuildCrosstalkMatrix(const std::array<TH2D*, 6>& xtalPairHists);

    void WriteCrosstalkMatrices(const std::string& fileName, const std::vector<TMatrixD>& xtalkMatrices);

    std::vector<TMatrixD> LoadCrosstalkMatrices(const std::string& fileName);

    std::vector<std::function<std::array<double, 4>(std::array<double, 4>)>> MakeCorrections(const std::string& xtalkCorrDir);

} // namespace CACrosstalkCorrection

#endif // CACROSSTALKCORRECTION_HPP