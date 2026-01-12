#include "qcovariate_loader.h"

#include <cmath>
#include <filesystem>
#include <fstream>
#include <ranges>
#include <vector>

#include <fmt/ranges.h>
#include <Eigen/Core>

#include "gelex/exception.h"
#include "gelex/logger.h"

// Internal
#include "../src/data/parser.h"

namespace gelex::detail
{
// =============================================================================
// QcovarLoader
// =============================================================================

QuantitativeCovariateLoader::QuantitativeCovariateLoader(
    const std::filesystem::path& path,
    bool iid_only)
{
    auto file = detail::open_file<std::ifstream>(path, std::ios::in);
    try
    {
        set_names(file);
        set_data(file, iid_only);
    }
    catch (const GelexException& e)
    {
        throw FileFormatException(
            std::format("{}:{}", path.string(), e.what()));
    }
    gelex::logging::get()->info(
        "Loaded {} samples with {} qcovars.", data_.size(), names_.size());
}

void QuantitativeCovariateLoader::set_names(std::ifstream& file)
{
    std::string line;
    std::getline(file, line);

    auto header = parse_header(line);
    if (header.size() < 3)
    {
        throw ColumnRangeException("Qcovar must have > 2 columns");
    }
    names_.clear();
    for (size_t i = 2; i < header.size(); ++i)
    {
        names_.emplace_back(header[i]);
    }
}

void QuantitativeCovariateLoader::set_data(std::ifstream& file, bool iid_only)
{
    std::string line;
    int n_line = 0;
    data_.reserve(1024);
    std::vector<double> values_buffer;

    while (std::getline(file, line))
    {
        n_line++;
        if (line.empty())
        {
            continue;
        }
        try
        {
            parse_all_doubles(line, values_buffer, 2);
            if (values_buffer.size() < names_.size())
            {
                throw DataParseException(
                    std::format(
                        "expected {} quantitative covariate values, but found "
                        "{}",
                        names_.size(),
                        values_buffer.size()));
            }
            // Check for nan/inf values
            bool has_invalid = false;
            for (double val : values_buffer)
            {
                if (std::isnan(val) || std::isinf(val))
                {
                    has_invalid = true;
                    break;
                }
            }
            if (has_invalid)
            {
                continue;
            }
            data_.emplace(parse_id(line, iid_only), values_buffer);
        }
        catch (const GelexException& e)
        {
            throw DataParseException(
                std::format("{}: {}", n_line + 1, e.what()));
        }
    }
}

QuantitativeCovariate QuantitativeCovariateLoader::load(
    const std::unordered_map<std::string, Eigen::Index>& id_map) const
{
    const auto n_samples = static_cast<Eigen::Index>(id_map.size());
    const auto n_covars = static_cast<Eigen::Index>(names_.size());

    Eigen::MatrixXd result(n_samples, n_covars);
    result.setConstant(std::numeric_limits<double>::quiet_NaN());

    for (const auto& [id, values] : data_)
    {
        if (auto it = id_map.find(id); it != id_map.end())
        {
            const Eigen::Index row_idx = it->second;
            result.row(row_idx)
                = Eigen::Map<const Eigen::RowVectorXd>(values.data(), n_covars);
        }
    }

    return QuantitativeCovariate{.names = names_, .X = std::move(result)};
}

}  // namespace gelex::detail
