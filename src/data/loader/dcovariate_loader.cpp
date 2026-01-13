#include "dcovariate_loader.h"

#include <algorithm>
#include <filesystem>
#include <fstream>
#include <map>
#include <string_view>
#include <vector>

#include <fmt/ranges.h>
#include <Eigen/Core>

#include "gelex/exception.h"

// Internal
#include "../src/data/parser.h"

namespace gelex::detail
{

namespace
{
auto iequals(std::string_view a, std::string_view b) -> bool
{
    return std::ranges::equal(
        a,
        b,
        [](unsigned char c1, unsigned char c2)
        { return std::tolower(c1) == std::tolower(c2); });
}

auto is_nan_or_inf_string(std::string_view sv) -> bool
{
    if (sv.size() < 3 || sv.size() > 4)
    {
        return false;
    }
    return iequals(sv, "nan") || iequals(sv, "inf") || iequals(sv, "+inf")
           || iequals(sv, "-inf");
}
}  // namespace

DiscreteCovariateLoader::DiscreteCovariateLoader(
    const std::filesystem::path& path,
    bool iid_only)
{
    auto file = open_file<std::ifstream>(path, std::ios::in);

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

auto DiscreteCovariateLoader::init_columns(std::ifstream& file) -> void
{
    std::string line;
    std::getline(file, line);

    auto header = parse_header(line);
    if (header.size() < 3)
    {
        throw ColumnRangeException(
            "categorical covariates must have > 2 columns");
    }
    column_names_.assign(header.begin() + 2, header.end());
    columns_.resize(column_names_.size());
}

auto DiscreteCovariateLoader::fill_columns(std::ifstream& file, bool iid_only)
    -> void
{
    std::string line;
    size_t n_line = 1;
    std::vector<std::string_view> buffer;
    buffer.reserve(column_names_.size());

    while (std::getline(file, line))
    {
        ++n_line;
        if (line.empty())
        {
            continue;
        }

        try
        {
            parse_string(line, buffer, 2);

            if (buffer.size() != column_names_.size())
            {
                throw DataParseException("Column count mismatch");
            }

            if (std::ranges::any_of(buffer, is_nan_or_inf_string))
            {
                continue;
            }

            sample_ids_.emplace_back(parse_id(line, iid_only));
            for (size_t i = 0; i < column_names_.size(); ++i)
            {
                columns_[i].data.push_back(
                    get_or_add_level(columns_[i], buffer[i]));
            }
        }
        catch (const GelexException& e)
        {
            throw DataParseException(std::format("{}: {}", n_line, e.what()));
        }
    }
}

auto DiscreteCovariateLoader::load(
    const std::unordered_map<std::string, Eigen::Index>& id_map) const
    -> DiscreteCovariate
{
    std::vector<Eigen::Index> file_indices;
    std::vector<Eigen::Index> target_indices;
    file_indices.reserve(id_map.size());
    target_indices.reserve(id_map.size());

    for (size_t i = 0; i < sample_ids_.size(); ++i)
    {
        if (auto it = id_map.find(sample_ids_[i]); it != id_map.end())
        {
            file_indices.push_back(static_cast<Eigen::Index>(i));
            target_indices.push_back(it->second);
        }
    }

    struct ColMeta
    {
        std::vector<int> global_id_to_local_rank;
        Eigen::Index start_col{};
        std::vector<std::string> sorted_levels;
    };
    std::vector<ColMeta> metas(columns_.size());
    Eigen::Index total_cols = 0;

    for (size_t i = 0; i < columns_.size(); ++i)
    {
        const auto& col = columns_[i];
        std::map<std::string, uint16_t> active_levels;

        for (auto row_idx : file_indices)
        {
            uint16_t val_id = col.data[row_idx];
            active_levels.emplace(col.levels[val_id], val_id);
        }

        if (active_levels.size() < 2)
        {
            continue;
        }

        metas[i].start_col = total_cols;
        metas[i].global_id_to_local_rank.assign(
            col.levels.size(), -2);  // -2: unused

        int rank = -1;  // -1: Reference, 0: First Dummy, 1: Second...
        for (const auto& [name, global_id] : active_levels)
        {
            metas[i].sorted_levels.push_back(name);
            metas[i].global_id_to_local_rank[global_id] = rank++;
        }
        total_cols += (static_cast<Eigen::Index>(active_levels.size()) - 1);
    }

    Eigen::MatrixXd X(id_map.size(), total_cols);
    X.setZero();

    std::vector<std::vector<std::string>> res_levels(columns_.size());
    std::vector<std::string> res_refs(columns_.size());

    for (size_t j = 0; j < columns_.size(); ++j)
    {
        if (metas[j].sorted_levels.empty())
        {
            continue;
        }
        res_refs[j] = metas[j].sorted_levels[0];
        res_levels[j] = std::move(metas[j].sorted_levels);

        const auto& mapping = metas[j].global_id_to_local_rank;
        const auto offset = metas[j].start_col;
        const auto& data_vec = columns_[j].data;

        for (size_t k = 0; k < file_indices.size(); ++k)
        {
            auto file_row = file_indices[k];
            auto target_row = target_indices[k];

            int rank = mapping[data_vec[file_row]];
            if (rank >= 0)
            {
                X(target_row, offset + rank) = 1.0;
            }
        }
    }

    return DiscreteCovariate{
        .names = column_names_,
        .levels = std::move(res_levels),
        .reference_levels = std::move(res_refs),
        .X = std::move(X)};
}

auto DiscreteCovariateLoader::get_or_add_level(
    ColumnData& column,
    std::string_view level) -> uint16_t
{
    if (auto it = column.level_map.find(level); it != column.level_map.end())
    {
        return it->second;
    }
    auto id = static_cast<uint16_t>(column.levels.size());
    column.levels.emplace_back(level);
    column.level_map.emplace(level, id);
    return id;
}

}  // namespace gelex::detail
