#include "covariate_processor.h"

#include <cmath>
#include <filesystem>
#include <format>
#include <fstream>
#include <sstream>
#include <string>

#include "../src/data/parser.h"

namespace gelex::detail
{

CovariateProcessor::CovariateProcessor(
    const std::filesystem::path& param_file_path)
{
    auto [intercept, continuous_coeffs, categorical_coeffs]
        = parse_param_file(param_file_path);

    intercept_ = intercept;
    continuous_coeffs_ = std::move(continuous_coeffs);
    categorical_coeffs_ = std::move(categorical_coeffs);
}

auto CovariateProcessor::parse_param_file(
    const std::filesystem::path& file_path)
    -> std::tuple<
        double,
        std::map<std::string, double>,
        std::map<std::string, std::map<std::string, double>>>
{
    auto file = detail::open_file<std::ifstream>(file_path, std::ios_base::in);

    std::string line;
    // Skip header line
    if (!std::getline(file, line))
    {
        throw InvalidFileException(
            std::format(
                "Parameter file '{}' is empty or has no header",
                file_path.string()));
    }

    double intercept = std::numeric_limits<double>::quiet_NaN();
    std::map<std::string, double> continuous_coeffs;
    std::map<std::string, std::map<std::string, double>> categorical_coeffs;

    // Process each line
    while (std::getline(file, line))
    {
        if (line.empty())
        {
            continue;
        }

        std::istringstream iss(line);
        std::string term;
        double mean{};
        double stddev{};
        double percentile_5{};
        double percentile_95{};
        double ess{};
        double rhat{};
        if (!(iss >> term >> mean >> stddev >> percentile_5 >> percentile_95
              >> ess >> rhat))
        {
            continue;  // Skip malformed lines
        }

        // Parse the term name and store coefficient
        parse_flat_name(
            term, mean, intercept, continuous_coeffs, categorical_coeffs);
    }

    // Validate that we found an intercept
    if (std::isnan(intercept))
    {
        throw InvalidDataException(
            std::format(
                "No intercept term found in parameter file '{}'",
                file_path.string()));
    }

    return std::make_tuple(
        intercept, std::move(continuous_coeffs), std::move(categorical_coeffs));
}

void CovariateProcessor::parse_flat_name(
    const std::string& flat_name,
    double coefficient,
    double& intercept,
    std::map<std::string, double>& continuous_coeffs,
    std::map<std::string, std::map<std::string, double>>& categorical_coeffs)
{
    if (flat_name == "Intercept")
    {
        intercept = coefficient;
        return;
    }

    // Check for categorical variable pattern (e.g., "Group_A", "Group_B")
    size_t underscore_pos = flat_name.find('_');
    if (underscore_pos != std::string::npos)
    {
        std::string var_name = flat_name.substr(0, underscore_pos);
        std::string category = flat_name.substr(underscore_pos + 1);
        categorical_coeffs[var_name][category] = coefficient;
    }
    else
    {
        // Treat as continuous variable
        continuous_coeffs[flat_name] = coefficient;
    }
}

double CovariateProcessor::predict(const IndividualData& data) const
{
    double score = intercept_;

    // Add continuous variable contributions
    for (const auto& [var_name, value] : data.continuous_values)
    {
        auto it = continuous_coeffs_.find(var_name);
        if (it != continuous_coeffs_.end())
        {
            score += value * it->second;
        }
    }

    // Add categorical variable contributions
    for (const auto& [var_name, category] : data.categorical_values)
    {
        auto var_it = categorical_coeffs_.find(var_name);
        if (var_it != categorical_coeffs_.end())
        {
            auto cat_it = var_it->second.find(category);
            if (cat_it != var_it->second.end())
            {
                score += cat_it->second;
            }
        }
    }

    return score;
}

}  // namespace gelex::detail
