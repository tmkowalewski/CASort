#ifndef CACROSSTALKCORRECTION_HPP
#define CACROSSTALKCORRECTION_HPP

// C++ Includes
#include <vector>

// ROOT Includes
#include <TGraphErrors.h>
#include <TH2D.h>


// Project Includes


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
        double alpha_xy = std::numeric_limits<double>::quiet_NaN();  // Crosstalk coefficient from x to y
        double alpha_yx = std::numeric_limits<double>::quiet_NaN();  // Crosstalk coefficient from y to x
        double alpha_xy_err = std::numeric_limits<double>::quiet_NaN();
        double alpha_yx_err = std::numeric_limits<double>::quiet_NaN();
        double chi2 = std::numeric_limits<double>::quiet_NaN();
        double ndf = std::numeric_limits<double>::quiet_NaN();
    };

    const double kTargetEnergy = 5018.98;       // Energy in keV of gamma ray used for crosstalk calibration
    const double kEnergyWindow = 15.0;          // Half-width of energy window in keV around target energy
    const double kEnergyCut = 244.0;            // Minimum energy in keV for cut in crosstalk plots
    const double kFitWindow = 30.0;             // Width of fit window in keV around target energy
    const unsigned int kMinCountsPerBin = 6;    // Minimum counts per x-bin slice for fitting

    // Function used to model crosstalk effect
    double CrosstalkFitFunction(double* x, double* par);

    std::shared_ptr<TGraphErrors> BuildCrosstalkGraph(const TH2D* hist);

    CrosstalkFit FitCrosstalkGraph(std::vector<TH2D*> xtal_pair_hists);

    std::vector<std::vector<std::function<double(double)>>> MakeCorrection();


} // namespace CACrosstalkCorrection

#endif // CACROSSTALKCORRECTION_HPP