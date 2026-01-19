#include "qcovariate_loader.h"

#include <cmath>
#include <filesystem>
#include <fstream>
#include <vector>

#include <fmt/ranges.h>
#include <Eigen/Core>

#include "gelex/exception.h"

// Internal
#include "../src/data/parser.h"

namespace gelex::detail
{

QuantitativeCovariateLoader::QuantitativeCovariateLoader(
    const std::filesystem::path& path,
    bool iid_only)
{
    auto file = detail::open_file<std::ifstream>(path, std::ios::in);
    try
    {
        init_columns(file);
        fill_columns(file, iid_only);
    }
    catch (const GelexException& e)
    {
        throw FileFormatException(
            std::format("{}:{}", path.string(), e.what()));
    }
}

auto QuantitativeCovariateLoader::init_columns(std::ifstream& file) -> void
{
    std::string line;
    std::getline(file, line);

    auto header = parse_header(line);
    if (header.size() <= IdColumnCount)
    {
        throw ColumnRangeException("Qcovar must have > 2 columns");
    }

    // Skip FID and IID columns (first 2 columns)
    column_names_.assign(header.begin() + IdColumnCount, header.end());
    columns_.resize(column_names_.size());
}

auto QuantitativeCovariateLoader::fill_columns(
    std::ifstream& file,
    bool iid_only) -> void
{
    std::string line;
    int line_number = 1;
    std::vector<double> values_buffer;
    const size_t n_covars = column_names_.size();

    while (std::getline(file, line))
    {
        line_number++;
        if (line.empty())
        {
            continue;
        }
        try
        {
            parse_all_doubles(line, values_buffer, IdColumnCount);
            if (values_buffer.size() != n_covars)
            {
                throw DataParseException(
                    std::format(
                        "expected {} quantitative covariate values, but found "
                        "{}",
                        n_covars,
                        values_buffer.size()));
            }

            // Skip rows with invalid values (NaN or Inf)
            if (std::ranges::any_of(
                    values_buffer,
                    [](double val) { return !is_valid_covariate_value(val); }))
            {
                continue;
            }

            sample_ids_.emplace_back(parse_id(line, iid_only));
            for (size_t i = 0; i < n_covars; ++i)
            {
                columns_[i].push_back(values_buffer[i]);
            }
        }
        catch (const GelexException& e)
        {
            throw DataParseException(
                std::format("{}: {}", line_number, e.what()));
        }
    }
}

auto QuantitativeCovariateLoader::load(
    const std::unordered_map<std::string, Eigen::Index>& id_map) const
    -> QuantitativeCovariate
{
    std::vector<Eigen::Index> file_indices;
    std::vector<Eigen::Index> target_indices;
    file_indices.reserve(id_map.size());
    target_indices.reserve(id_map.size());

    // Build index mapping from sample file IDs to target position
    for (size_t i = 0; i < sample_ids_.size(); ++i)
    {
        if (auto it = id_map.find(sample_ids_[i]); it != id_map.end())
        {
            file_indices.push_back(static_cast<Eigen::Index>(i));
            target_indices.push_back(it->second);
        }
    }

    const auto n_samples = static_cast<Eigen::Index>(id_map.size());
    const auto n_covars = static_cast<Eigen::Index>(column_names_.size());

    Eigen::MatrixXd X(n_samples, n_covars);
    X.setConstant(std::numeric_limits<double>::quiet_NaN());

    // Fill matrix with covariate values in row-major access pattern
    // for better cache locality
    for (size_t match_idx = 0; match_idx < file_indices.size(); ++match_idx)
    {
        const auto file_row = file_indices[match_idx];
        const auto target_row = target_indices[match_idx];

        for (size_t covar_idx = 0; covar_idx < columns_.size(); ++covar_idx)
        {
            X(target_row, static_cast<Eigen::Index>(covar_idx))
                = columns_[covar_idx][file_row];
        }
    }

    return QuantitativeCovariate{.names = column_names_, .X = std::move(X)};
}

}  // namespace gelex::detail

auto gelex::detail::QuantitativeCovariateLoader::is_valid_covariate_value(
    double value) -> bool
{
    return !std::isnan(value) && !std::isinf(value);
}
