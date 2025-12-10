#include "covar_effect_loader.h"

#include "gelex/exception.h"
#include "gelex/logger.h"

#include <cmath>
#include <filesystem>
#include <format>
#include <fstream>
#include <sstream>
#include <string>

#include "../src/data/parser.h"

namespace gelex::detail
{

CovarEffectLoader::CovarEffectLoader(
    const std::filesystem::path& param_file_path)
{
    auto [intercept, continuous_coeffs, categorical_coeffs]
        = parse_param_file(param_file_path);

    intercept_ = intercept;
    continuous_coeffs_ = std::move(continuous_coeffs);
    categorical_coeffs_ = std::move(categorical_coeffs);

    auto logger = gelex::logging::get();
    size_t categorical_categories = 0;
    for (const auto& [var, categories] : categorical_coeffs_)
    {
        categorical_categories += categories.size();
    }
    logger->info(
        "Loaded covariate effects: intercept={}, continuous vars={}, "
        "categorical vars={} categories",
        intercept_,
        continuous_coeffs_.size(),
        categorical_categories);
}

auto CovarEffectLoader::parse_param_file(const std::filesystem::path& file_path)
    -> std::tuple<
        double,
        std::map<std::string, double>,
        std::map<std::string, std::map<std::string, double>>>
{
    auto file = detail::open_file<std::ifstream>(file_path, std::ios::in);

    std::string line;
    // Skip header line
    if (!std::getline(file, line))
    {
        throw FileFormatException(
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
        throw DataParseException(
            std::format(
                "No intercept term found in parameter file '{}'",
                file_path.string()));
    }

    return std::make_tuple(
        intercept, std::move(continuous_coeffs), std::move(categorical_coeffs));
}

void CovarEffectLoader::parse_flat_name(
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

double CovarEffectLoader::predict(const IndividualData& data) const
{
    double score = intercept_;

    // Add continuous variable contributions
    for (const auto& [var_name, value] : data.continuous_values)
    {
        auto it = continuous_coeffs_.find(var_name);
        if (it == continuous_coeffs_.end())
        {
            throw InvalidInputException(
                std::format("unknown continuous variable '{}'", var_name));
        }
        score += value * it->second;
    }

    // Add categorical variable contributions
    for (const auto& [var_name, category] : data.categorical_values)
    {
        auto var_it = categorical_coeffs_.find(var_name);
        if (var_it == categorical_coeffs_.end())
        {
            throw InvalidInputException(
                std::format("unknown categorical variable '{}'", var_name));
        }

        auto cat_it = var_it->second.find(category);
        if (cat_it == var_it->second.end())
        {
            throw InvalidInputException(
                std::format(
                    "unknown category '{}' for categorical variable '{}'",
                    category,
                    var_name));
        }
        score += cat_it->second;
    }

    return score;
}

}  // namespace gelex::detail
