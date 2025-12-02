#include "ccovariate_loader.h"

#include <algorithm>
#include <filesystem>
#include <fstream>
#include <ranges>
#include <string_view>
#include <unordered_set>
#include <vector>

#include <fmt/ranges.h>
#include <Eigen/Core>

#include "gelex/exception.h"
#include "gelex/logger.h"

// Internal
#include "../src/data/parser.h"

namespace gelex::detail
{

CCovarLoader::CCovarLoader(const std::filesystem::path& path, bool iid_only)
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
        "Loaded {} samples with {} categorical covars.",
        raw_data_.size(),
        names_.size());
}

void CCovarLoader::set_names(std::ifstream& file)
{
    std::string line;
    std::getline(file, line);

    auto header = parse_header(line);
    if (header.size() < 3)
    {
        throw ColumnRangeException(
            "categorical covariates must have > 2 columns");
    }
    names_.clear();
    for (size_t i = 2; i < header.size(); ++i)
    {
        names_.emplace_back(header[i]);
    }
}

void CCovarLoader::set_data(std::ifstream& file, bool iid_only)
{
    std::string line;
    int n_line = 0;
    std::vector<std::string_view> value_buffer;
    value_buffer.reserve(names_.size());

    while (std::getline(file, line))
    {
        n_line++;
        if (line.empty())
        {
            continue;
        }

        try
        {
            parse_string(line, value_buffer, 2);

            if (value_buffer.size() != names_.size())
            {
                throw DataParseException("Column count mismatch");
            }

            raw_data_.emplace(
                parse_id(line, iid_only),
                std::ranges::to<std::vector<std::string>>(value_buffer));
        }
        catch (const GelexException& e)
        {
            throw DataParseException(
                std::format("{}: {}", n_line + 1, e.what()));
        }
    }
}

Eigen::MatrixXd CCovarLoader::load(
    const std::unordered_map<std::string, Eigen::Index>& id_map) const
{
    auto [valid_ids, levels_per_col] = get_valid_samples_and_levels(id_map);
    auto [local_encodings, total_cols] = build_local_encodings(levels_per_col);

    return fill_matrix(id_map, valid_ids, local_encodings, total_cols);
}

CCovarLoader::IntersectResult CCovarLoader::get_valid_samples_and_levels(
    const std::unordered_map<std::string, Eigen::Index>& id_map) const
{
    std::vector<std::string_view> valid_ids;
    valid_ids.reserve(id_map.size());

    std::vector<std::unordered_set<std::string_view>> levels_per_col(
        names_.size());

    for (const auto& [id, row_values] : raw_data_)
    {
        if (id_map.contains(id))
        {
            valid_ids.emplace_back(id);

            for (size_t i = 0;
                 i < row_values.size() && i < levels_per_col.size();
                 ++i)
            {
                if (!row_values[i].empty())
                {
                    levels_per_col[i].insert(row_values[i]);
                }
            }
        }
    }
    return {std::move(valid_ids), std::move(levels_per_col)};
}

CCovarLoader::EncodingResult CCovarLoader::build_local_encodings(
    const std::vector<std::unordered_set<std::string_view>>& levels_per_col)
    const
{
    std::vector<std::unordered_map<std::string_view, EncodedCovariate>>
        encodings(names_.size());
    Eigen::Index total_dummy_vars = 0;

    for (size_t i = 0; i < names_.size(); ++i)
    {
        const auto& levels = levels_per_col[i];

        if (levels.size() < 2)
        {
            continue;
        }

        std::vector<std::string_view> sorted_levels(
            levels.begin(), levels.end());
        std::ranges::sort(sorted_levels);

        Eigen::Index current_offset = (total_dummy_vars);

        encodings[i].emplace(
            sorted_levels[0], EncodedCovariate{-1, current_offset});

        for (Eigen::Index j = 1;
             j < static_cast<Eigen::Index>(sorted_levels.size());
             ++j)
        {
            encodings[i].emplace(
                sorted_levels[j], EncodedCovariate{j - 1, current_offset});
        }

        total_dummy_vars
            += (static_cast<Eigen::Index>(sorted_levels.size()) - 1);
    }

    return {std::move(encodings), total_dummy_vars};
}

Eigen::MatrixXd CCovarLoader::fill_matrix(
    const std::unordered_map<std::string, Eigen::Index>& id_map,
    const std::vector<std::string_view>& valid_ids,
    const std::vector<std::unordered_map<std::string_view, EncodedCovariate>>&
        local_encodings,
    Eigen::Index total_cols) const
{
    Eigen::MatrixXd result(id_map.size(), total_cols);
    result.setZero();

    for (const auto& id_sv : valid_ids)
    {
        std::string id(id_sv);

        const Eigen::Index row_idx = id_map.at(id);
        const auto& raw_values = raw_data_.at(id);

        const size_t n_cols_process
            = std::min(raw_values.size(), local_encodings.size());

        for (size_t i = 0; i < n_cols_process; ++i)
        {
            const auto& encoding_map = local_encodings[i];

            if (encoding_map.empty())
            {
                continue;
            }

            auto enc_it = encoding_map.find(raw_values[i]);
            if (enc_it != encoding_map.end())
            {
                const auto& info = enc_it->second;

                if (info.active_index >= 0)
                {
                    result(row_idx, info.start_col_offset + info.active_index)
                        = 1.0;
                }
            }
        }
    }

    return result;
}

}  // namespace gelex::detail
