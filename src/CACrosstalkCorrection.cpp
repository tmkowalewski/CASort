// C++ Includes

// ROOT Includes
#include <TH2D.h>
#include <TF1.h>

// Project Includes
#include "CACrosstalkCorrection.hpp"

double CACrosstalkCorrection::CrosstalkFitFunction(double* x, double* par)
{
    const double alpha_xy = par[0];
    const double alpha_yx = par[1];
    const double E_gamma = par[2];

    const double k0 = alpha_xy / (1.0 - alpha_xy);
    const double k1 = alpha_yx / (1.0 - alpha_yx);

    const double slope = -(1.0 + k0) / (1.0 + k1);
    const double intercept = E_gamma * (1.0 + k0 / (1.0 + k1));
    return x[0] * slope + intercept;
}

std::shared_ptr<TGraphErrors> CACrosstalkCorrection::BuildCrosstalkGraph(const TH2D* hist)
{
    // Get details from histogram
    const size_t x_bins = hist->GetXaxis()->GetNbins();
    const size_t y_bins = hist->GetYaxis()->GetNbins();
    const auto hist_xax = hist->GetXaxis();
    const auto hist_yax = hist->GetYaxis();

    // Make a new plot for the crosstalk graph
    auto graph = std::make_shared<TGraphErrors>();
    try
    {
        // Set the title to EivEj_xtalk
        graph->SetNameTitle(Form("%sv%s_xtalk", hist_xax->GetName(), hist_yax->GetName()),
            Form("%sv%s_xtalk;E%d Measured Energy (keV); E%d Measured Energy (keV)",
                hist_xax->GetName(), hist_yax->GetName(),
                hist_xax->GetTitle()[std::string(hist_xax->GetTitle()).find('E') + 1],
                hist_yax->GetTitle()[std::string(hist_yax->GetTitle()).find('E') + 1]));
    }
    catch (...)
    {
        graph->SetNameTitle("xtalk_graph", "Crosstalk Graph;E_x Measured Energy (keV); E_y Measured Energy (keV)");
    }

    for (size_t ix = 1; ix <= x_bins; ++ix)
    {
        const double energy_x = hist_xax->GetBinCenter(ix);
        const double y_low = kTargetEnergy - kEnergyWindow / 2.0 - energy_x;  // y = E - x for m = 2 add-back events
        const double y_high = kTargetEnergy + kEnergyWindow / 2.0 - energy_x;

        size_t iy_min = std::max(1, hist_yax->FindBin(y_low));
        size_t iy_max = std::min(static_cast<int>(y_bins), hist_yax->FindBin(y_high));

        double sum_weights, sum_weighted_energy_y, sum_weighted_energy_y2;

        for (size_t iy = iy_min; iy <= iy_max; ++iy)
        {
            const double bin_content = hist->GetBinContent(ix, iy);
            if (bin_content < kMinCountsPerBin)
                continue;

            const double energy_y = hist_yax->GetBinCenter(iy);
            sum_weights += bin_content;
            sum_weighted_energy_y += bin_content * energy_y;
            sum_weighted_energy_y2 += bin_content * energy_y * energy_y;
        }

        const auto mean_Ey = sum_weighted_energy_y / sum_weights;
        const auto var_Ey = std::max(0.0, (sum_weighted_energy_y2 / sum_weights) - (mean_Ey * mean_Ey));
        const auto err_Ey = std::sqrt(var_Ey / sum_weights);

        const auto n_points = graph->GetN();
        graph->SetPoint(n_points, energy_x, mean_Ey);
        graph->SetPointError(n_points, hist_xax->GetBinWidth(ix) / 2.0, err_Ey);
    }

    return graph;
}

CACrosstalkCorrection::CrosstalkFit CACrosstalkCorrection::FitCrosstalkGraph(const TH2D* const hist)
{
    auto graph = BuildCrosstalkGraph(hist);

    auto fit_func = std::make_unique<TF1>("crosstalk_fit_func", CrosstalkFitFunction, kTargetEnergy - kFitWindow / 2.0, kTargetEnergy + kFitWindow / 2.0, 3);

    fit_func->SetParNames("alpha_xy", "alpha_yx", "E_gamma");
    fit_func->SetParameters(1e-4, 1e-4, kTargetEnergy);
    fit_func->FixParameter(2, kTargetEnergy);  // Fix E_gamma to known value

    graph->Fit(fit_func.get(), "QRS"); // Q quiet, R range, S store

    CrosstalkFit result;
    result.valid = true;
    result.alpha_xy = fit_func->GetParameter(0);
    result.alpha_yx = fit_func->GetParameter(1);
    result.alpha_xy_err = fit_func->GetParError(0);
    result.alpha_yx_err = fit_func->GetParError(1);
    result.chi2 = fit_func->GetChisquare();
    result.ndf = fit_func->GetNDF();

    return result;
}

std::vector<std::vector<std::function<double(double)>>> CACrosstalkCorrection::MakeCorrection()
{

}



