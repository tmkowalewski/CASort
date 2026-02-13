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
    const double alphaXY = par[0];
    const double alphaYX = par[1];
    const double gammaE = par[2];

    const double k0 = alphaXY / (1.0 - alphaXY);
    const double k1 = alphaYX / (1.0 - alphaYX);

    const double slope = -(1.0 + k0) / (1.0 + k1);
    const double intercept = gammaE * (1.0 + k0 / (1.0 + k1));
    return x[0] * slope + intercept;
}

void CACrosstalkCorrection::FillXTalkHistograms(const std::array<std::shared_ptr<TH2D>, 6>& xtalkPairHists, const std::array<double, 4>& xtalE, std::array<double, 4>& xtalT)
{
    static std::array<std::pair<short, short>, 6> xtalPairs = {{{0, 1}, {0, 2}, {0, 3}, {1, 2}, {1, 3}, {2, 3}}};
    for (size_t i = 0; i < xtalPairs.size(); i++)
    {
        auto [xtalX, xtalY] = xtalPairs[i];
        if (!std::isnan(xtalE[xtalX]) && !std::isnan(xtalE[xtalY]) && (fabs(xtalT[xtalX] - xtalT[xtalY]) < CAAddBack::kAddBackWindow))
        {
            xtalkPairHists[i]->Fill(xtalE[xtalX], xtalE[xtalY]);
        }
    }
}

std::shared_ptr<TGraphErrors> CACrosstalkCorrection::BuildCrosstalkGraph(const TH2D* hist)
{
    // Get details from histogram
    const size_t nBinsX = hist->GetXaxis()->GetNbins();
    const size_t nBinsY = hist->GetYaxis()->GetNbins();
    const auto histXaxis = hist->GetXaxis();
    const auto histYaxis = hist->GetYaxis();

    // Make a new plot for the crosstalk graph
    auto graph = std::make_shared<TGraphErrors>();
    try
    {
        // Set the title to EivEj_xtalk
        graph->SetNameTitle(Form("%s_gr", hist->GetName()), Form("%s;%s;%s", hist->GetTitle(), histXaxis->GetTitle(), histYaxis->GetTitle()));
    }
    catch (...)
    {
        graph->SetNameTitle("xtalk_graph", "Crosstalk Graph;E_x Measured Energy (keV); E_y Measured Energy (keV)");
    }

    for (size_t ix = 1; ix <= nBinsX; ++ix)
    {
        const double energyX = histXaxis->GetBinCenter(ix);
        const double yLow = std::max(kTargetEnergy - kEnergyWindow / 2.0 - energyX, 0.0);
        const double yHigh = std::max(kTargetEnergy + kEnergyWindow / 2.0 - energyX, 0.0);

        size_t iyMin = std::max(1, histYaxis->FindBin(yLow));
        size_t iyMax = std::min(static_cast<int>(nBinsY), histYaxis->FindBin(yHigh));

        double sumWeights = 0.0, sumWeightedEnergyY = 0.0, sumWeightedEnergyY2 = 0.0;

        for (size_t iy = iyMin; iy <= iyMax; ++iy)
        {
            const double binContent = hist->GetBinContent(ix, iy);
            if (binContent < kMinCountsPerBin)
                continue;
            const double energyY = histYaxis->GetBinCenter(iy);
            sumWeights += binContent;
            sumWeightedEnergyY += binContent * energyY;
            sumWeightedEnergyY2 += binContent * energyY * energyY;
        }

        if (sumWeights < kMinCountsPerBin)
            continue;
        const auto meanEY = sumWeightedEnergyY / sumWeights;
        const auto varEY = std::max(0.0, (sumWeightedEnergyY2 / sumWeights) - (meanEY * meanEY));
        const auto errEY = std::sqrt(varEY / sumWeights);

        // printf("X-bin %zu (Energy = %.1f keV, y-bin range = [%.1f, %.1f])\n", ix, energyX, histYaxis->GetBinCenter(iyMin), histYaxis->GetBinCenter(iyMax));
        // printf("Mean E_y = %.3f keV, Std Dev E_y = %.3f keV, Counts = %.3f\n", meanEY, std::sqrt(varEY), sumWeights);

        const auto nPoints = graph->GetN();
        graph->SetPoint(nPoints, energyX, meanEY);
        graph->SetPointError(nPoints, histXaxis->GetBinWidth(ix) / 2.0, errEY);
    }

    return graph;
}

CACrosstalkCorrection::CrosstalkFit CACrosstalkCorrection::FitCrosstalkCorrection(const TH2D* hist)
{
    auto graph = BuildCrosstalkGraph(hist);

    auto fitFunc = std::make_unique<TF1>("crosstalk_fit_func", CrosstalkFitFunction, 0, kTargetEnergy + kFitWindow / 2.0, 3);

    fitFunc->SetParNames("alphaXY", "alphaYX", "gammaE");
    fitFunc->SetParameters(1e-4, 1e-4, kTargetEnergy);
    fitFunc->FixParameter(2, kTargetEnergy); // Fix gammaE to known value

    graph->Fit(fitFunc.get(), "RS"); // Q quiet, R range, S store

    CrosstalkFit result;
    result.valid = true;
    result.alphaXY = fitFunc->GetParameter(0);
    result.alphaYX = fitFunc->GetParameter(1);
    result.alphaXYErr = fitFunc->GetParError(0);
    result.alphaYXErr = fitFunc->GetParError(1);
    result.chi2 = fitFunc->GetChisquare();
    result.ndf = fitFunc->GetNDF();

    return result;
}

TMatrixD CACrosstalkCorrection::BuildCrosstalkMatrix(const std::array<TH2D*, 6>& xtalPairHists)
{
    TMatrixD xtalkMatrix(4, 4);
    xtalkMatrix.Zero();

    auto xtalPairs = std::array<std::pair<short, short>, 6>{{{0, 1}, {0, 2}, {0, 3}, {1, 2}, {1, 3}, {2, 3}}};

    for (size_t i = 0; i < xtalPairHists.size(); ++i)
    {
        auto fitResult = FitCrosstalkCorrection(xtalPairHists[i]);
        if (!fitResult.valid)
        {
            throw std::runtime_error("[ERROR] Crosstalk fit failed for histogram: " + std::string(xtalPairHists[i]->GetName()));
        }

        auto [xtalX, xtalY] = xtalPairs[i];
        xtalkMatrix(xtalX, xtalY) = fitResult.alphaXY;
        xtalkMatrix(xtalY, xtalX) = fitResult.alphaYX;
    }

    return xtalkMatrix;
}

void CACrosstalkCorrection::WriteCrosstalkMatrices(const std::string& fileName, const std::vector<TMatrixD>& xtalkMatrices)
{
    FILE* outputFile = fopen(fileName.c_str(), "w");
    if (!outputFile)
    {
        throw std::runtime_error("[ERROR] Failed to open file for writing: " + fileName);
    }

    fprintf(outputFile, "# Channel\t a_i0\t a_i1\t a_i2\t a_i3\n");
    for (size_t i = 0; i < xtalkMatrices.size(); ++i)
    {
        fprintf(outputFile, "# Detector %zu\n", i);
        const auto& xtalkMatrix = xtalkMatrices[i];
        for (size_t i = 0; i < 4; ++i)
        {
            fprintf(outputFile, "%zu\t%14.10f\t%14.10f\t%14.10f\t%14.10f\n", i, xtalkMatrix(i, 0), xtalkMatrix(i, 1), xtalkMatrix(i, 2), xtalkMatrix(i, 3));
        }
    }
    fclose(outputFile);
}

std::vector<TMatrixD> CACrosstalkCorrection::LoadCrosstalkMatrices(const std::string& fileName)
{
    std::vector<TMatrixD> xtalkMatrices;

    FILE* inputFile = fopen(fileName.c_str(), "r");
    if (!inputFile)
    {
        throw std::runtime_error("[ERROR] Failed to open file for reading: " + fileName);
    }

    auto raw_data = CAUtilities::ReadCAFile(fileName);
    auto xtalkMatrix = TMatrixD(4, 4);
    xtalkMatrix.Zero();
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
                xtalkMatrix(channel, j) = row[j + 1];
            }
        }

        xtalkMatrices.push_back(xtalkMatrix);
    }

    return xtalkMatrices;
}

std::vector<std::function<std::array<double, 4>(std::array<double, 4>)>> CACrosstalkCorrection::MakeCorrections(const std::string& fileName)
{
    std::vector<std::function<std::array<double, 4>(std::array<double, 4>)>> corrections;

    auto xtalkMatrices = LoadCrosstalkMatrices(fileName);

    for (const auto& xtalkMatrix : xtalkMatrices)
    {
        corrections.push_back([xtalkMatrix](std::array<double, 4> measE) -> std::array<double, 4>
                              {
            std::array<double, 4> corrE;
            for (size_t i = 0; i < 4; ++i)
            {
                double correction = 0.0;
                for (size_t j = 0; j < 4; ++j)
                {
                    correction += xtalkMatrix(i, j) * measE[j];
                }
                corrE[i] = measE[i] - correction;
            }
            return corrE; });
    }

    return corrections;
}
