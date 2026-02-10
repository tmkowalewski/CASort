// C++ Includes
#include <fstream>

// ROOT Includes
#include <TF1.h>

// Project Includes
#include "CAAddBack.hpp"
#include "CACrosstalkCorrection.hpp"
#include "CAUtilities.hpp"

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

void CACrosstalkCorrection::FillXTalkHistograms(const std::array<std::shared_ptr<TH2D>, 6>& xtalk_pair_hists, const std::array<double, 4>& xtal_E, std::array<double, 4>& xtal_T)
{
    static std::array<std::pair<short, short>, 6> xtal_pairs = {{{0, 1}, {0, 2}, {0, 3}, {1, 2}, {1, 3}, {2, 3}}};
    for (size_t i = 0; i < xtal_pairs.size(); i++)
    {
        auto [xtal1, xtal2] = xtal_pairs[i];
        if (!std::isnan(xtal_E[xtal1]) && !std::isnan(xtal_E[xtal2]) && (fabs(xtal_T[xtal1] - xtal_T[xtal2]) < CAAddBack::kAddBackWindow))
        {
            xtalk_pair_hists[i]->Fill(xtal_E[xtal1], xtal_E[xtal2]);
        }
    }
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
        graph->SetNameTitle(Form("%s_gr", hist->GetName()), Form("%s;%s;%s", hist->GetTitle(), hist_xax->GetTitle(), hist_yax->GetTitle()));
    }
    catch (...)
    {
        graph->SetNameTitle("xtalk_graph", "Crosstalk Graph;E_x Measured Energy (keV); E_y Measured Energy (keV)");
    }

    for (size_t ix = 1; ix <= x_bins; ++ix)
    {
        const double energy_x = hist_xax->GetBinCenter(ix);
        const double y_low = std::max(kTargetEnergy - kEnergyWindow / 2.0 - energy_x, 0.0);
        const double y_high = std::max(kTargetEnergy + kEnergyWindow / 2.0 - energy_x, 0.0);

        size_t iy_min = std::max(1, hist_yax->FindBin(y_low));
        size_t iy_max = std::min(static_cast<int>(y_bins), hist_yax->FindBin(y_high));

        double sum_weights = 0.0, sum_weighted_energy_y = 0.0, sum_weighted_energy_y2 = 0.0;

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

        if (sum_weights < kMinCountsPerBin)
            continue;
        const auto mean_Ey = sum_weighted_energy_y / sum_weights;
        const auto var_Ey = std::max(0.0, (sum_weighted_energy_y2 / sum_weights) - (mean_Ey * mean_Ey));
        const auto err_Ey = std::sqrt(var_Ey / sum_weights);

        // printf("X-bin %zu (Energy = %.1f keV, y-bin range = [%.1f, %.1f])\n", ix, energy_x, hist_yax->GetBinCenter(iy_min), hist_yax->GetBinCenter(iy_max));
        // printf("Mean E_y = %.3f keV, Std Dev E_y = %.3f keV, Counts = %.3f\n", mean_Ey, std::sqrt(var_Ey), sum_weights);

        const auto n_points = graph->GetN();
        graph->SetPoint(n_points, energy_x, mean_Ey);
        graph->SetPointError(n_points, hist_xax->GetBinWidth(ix) / 2.0, err_Ey);
    }

    return graph;
}

CACrosstalkCorrection::CrosstalkFit CACrosstalkCorrection::FitCrosstalkCorrection(const TH2D* hist)
{
    auto graph = BuildCrosstalkGraph(hist);

    auto fit_func = std::make_unique<TF1>("crosstalk_fit_func", CrosstalkFitFunction, 0, kTargetEnergy + kFitWindow / 2.0, 3);

    fit_func->SetParNames("alpha_xy", "alpha_yx", "E_gamma");
    fit_func->SetParameters(1e-4, 1e-4, kTargetEnergy);
    fit_func->FixParameter(2, kTargetEnergy); // Fix E_gamma to known value

    graph->Fit(fit_func.get(), "RS"); // Q quiet, R range, S store

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

TMatrixD CACrosstalkCorrection::BuildCrosstalkMatrix(std::array<TH2D*, 6>& xtal_pair_hists)
{
    TMatrixD A(4, 4);
    A.Zero();

    auto xtal_pairs = std::array<std::pair<short, short>, 6>{{{0, 1}, {0, 2}, {0, 3}, {1, 2}, {1, 3}, {2, 3}}};

    for (size_t i = 0; i < xtal_pair_hists.size(); ++i)
    {
        auto fit_result = FitCrosstalkCorrection(xtal_pair_hists[i]);
        if (!fit_result.valid)
        {
            throw std::runtime_error("[ERROR] Crosstalk fit failed for histogram: " + std::string(xtal_pair_hists[i]->GetName()));
        }

        auto [xtal1, xtal2] = xtal_pairs[i];
        A(xtal1, xtal2) = fit_result.alpha_xy;
        A(xtal2, xtal1) = fit_result.alpha_yx;
    }

    return A;
}

void CACrosstalkCorrection::WriteCrosstalkMatrices(std::string filename, const std::vector<TMatrixD>& xtalk_matrices)
{
    FILE* out_file = fopen(filename.c_str(), "w");
    if (!out_file)
    {
        throw std::runtime_error("[ERROR] Failed to open file for writing: " + filename);
    }

    fprintf(out_file, "# Channel\t a_i0\t a_i1\t a_i2\t a_i3\n");
    for (size_t i = 0; i < xtalk_matrices.size(); ++i)
    {
        fprintf(out_file, "# Detector %zu\n", i);
        const auto& xtalk_matrix = xtalk_matrices[i];
        for (size_t i = 0; i < 4; ++i)
        {
            fprintf(out_file, "%zu\t%14.10f\t%14.10f\t%14.10f\t%14.10f\n", i, xtalk_matrix(i, 0), xtalk_matrix(i, 1), xtalk_matrix(i, 2), xtalk_matrix(i, 3));
        }
    }
    fclose(out_file);
}

std::vector<TMatrixD> CACrosstalkCorrection::LoadCrosstalkMatrices(const std::string& filename)
{
    std::vector<TMatrixD> matrices;

    FILE* in_file = fopen(filename.c_str(), "r");
    if (!in_file)
    {
        throw std::runtime_error("[ERROR] Failed to open file for reading: " + filename);
    }

    auto raw_data = CAUtilities::ReadCAFile(filename);
    auto M = TMatrixD(4, 4);
    M.Zero();
    for (const auto& matrix_data : raw_data)
    {
        for (const auto& row : matrix_data)
        {
            if (row.size() != 5)
            {
                throw std::runtime_error("[ERROR] Invalid row size in crosstalk matrix file. Expected 5 columns (channel, a_i0, a_i1, a_i2, a_i3)");
            }
            size_t channel = static_cast<size_t>(row[0]);
            if (channel >= 4)
            {
                throw std::runtime_error("[ERROR] Invalid column size in crosstalk matrix file. Expected columns of size 4 (a_0j, a_1j, a_2j, a_3j)");
            }
            for (size_t j = 0; j < 4; ++j)
            {
                M(channel, j) = row[j + 1];
            }
        }

        matrices.push_back(M);
    }

    return matrices;
}

std::vector<std::function<std::array<double, 4>(std::array<double, 4>)>> CACrosstalkCorrection::MakeCorrections(const std::string& filename)
{
    std::vector<std::function<std::array<double, 4>(std::array<double, 4>)>> corrections;

    auto xtalk_matrices = LoadCrosstalkMatrices(filename);

    for (const auto& M : xtalk_matrices)
    {
        corrections.push_back([M](std::array<double, 4> E_meas) -> std::array<double, 4>
                              {
            std::array<double, 4> E_corr;
            for (size_t i = 0; i < 4; ++i)
            {
                double correction = 0.0;
                for (size_t j = 0; j < 4; ++j)
                {
                    correction += M(i, j) * E_meas[j];
                }
                E_corr[i] = E_meas[i] - correction;
            }
            return E_corr; });
    }

    return corrections;
}
