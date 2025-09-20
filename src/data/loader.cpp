#include "loader.h"

#include <expected>
#include <filesystem>
#include <fstream>
#include <ranges>
#include <string>
#include <string_view>
#include <unordered_map>
#include <unordered_set>
#include <vector>

#include <fmt/ranges.h>
#include <Eigen/Core>

#include "../src/data/parser.h"
#include "gelex/error.h"
#include "gelex/logger.h"

namespace gelex::detail
{

auto PhenotypeLoader::create(
    std::string_view path,
    int pheno_column,
    bool iid_only) -> std::expected<PhenotypeLoader, Error>
{
    auto file = detail::openfile<std::ifstream>(path);
    if (!file)
    {
        return std::unexpected(file.error());
    }
    auto header = get_header(*file, pheno_column);
    if (!header)
    {
        return std::unexpected(
            enrich_with_file_info(std::move(header.error()), path));
    }
    size_t cols = header->size();
    auto map = read(*file, pheno_column, cols, iid_only);
    if (!map)
    {
        return std::unexpected(
            enrich_with_file_info(std::move(map.error()), path));
    }

    auto logger = gelex::logging::get();
    logger->info(
        "Loaded {} samples with phenotype '{}'",
        map->size(),
        (*header)[pheno_column]);

    return PhenotypeLoader(std::move((*header)[pheno_column]), std::move(*map));
}

void PhenotypeLoader::intersect(std::unordered_set<std::string>& id_set) const
{
    std::erase_if(
        id_set,
        [this](const std::string& current_id)
        { return !phenotype_data_.contains(current_id); });
}

Eigen::VectorXd PhenotypeLoader::load(
    const std::unordered_map<std::string, Eigen::Index>& id_map) const
{
    Eigen::VectorXd data(id_map.size(), 1);

    for (const auto& [id, value] : phenotype_data_)
    {
        if (auto it = id_map.find(id); it != id_map.end())
        {
            data(it->second) = value;
        }
    }
    return data;
}

auto PhenotypeLoader::get_header(std::ifstream& file, size_t target_column)
    -> std::expected<std::vector<std::string>, Error>
{
    std::string line;
    std::getline(file, line);
    const auto header_view = parse_header(line);

    if (!header_view)
    {
        return std::unexpected(header_view.error());
    }

    if (target_column < 2 || target_column >= (header_view->size()))
    {
        return std::unexpected(
            Error{
                ErrorCode::InvalidRange,
                std::format(
                    "Phenotype column index {} is out of range [2, {}) ",
                    target_column,
                    header_view->size())});
    }
    return std::ranges::to<std::vector<std::string>>(*header_view);
}

auto PhenotypeLoader::read(
    std::ifstream& file,
    size_t target_column,
    size_t expected_columns,
    bool iid_only)
    -> std::expected<std::unordered_map<std::string, double>, Error>
{
    std::unordered_map<std::string, double> phenotype_data;

    std::string line;
    int n_line{};
    while (std::getline(file, line))
    {
        if (line.empty())
        {
            continue;
        }

        if (count_num_columns(line) != expected_columns)
        {
            return std::unexpected(
                Error{
                    ErrorCode::InconsistColumnCount,
                    std::format(
                        "Inconsistent number of columns (line {})",
                        n_line + 2)});
        }

        if (auto id_str = parse_id(line, iid_only); id_str)
        {
            if (auto value = parse_nth_double(line, target_column))
            {
                phenotype_data.emplace(std::move(*id_str), *value);
            }
            else
            {
                return std::unexpected(enrich_with_line_info(
                    std::move(value.error()), n_line + 2));
            }
        }
        else
        {
            return std::unexpected(
                enrich_with_line_info(std::move(id_str.error()), n_line + 2));
        }
        n_line++;
    }
    return phenotype_data;
}

auto QcovarLoader::create(std::string_view path, bool iid_only)
    -> std::expected<QcovarLoader, Error>
{
    auto file = detail::openfile<std::ifstream>(path);
    if (!file)
    {
        return std::unexpected(file.error());
    }
    auto header = get_header(*file);
    if (!header)
    {
        return std::unexpected(
            enrich_with_file_info(std::move(header.error()), path));
    }
    size_t cols = header->size();
    auto map = read(*file, cols, iid_only);
    if (!map)
    {
        return std::unexpected(
            enrich_with_file_info(std::move(map.error()), path));
    }

    std::vector<std::string> covar_names;
    covar_names.reserve(cols - 2);
    for (size_t i = 2; i < cols; ++i)
    {
        covar_names.push_back(std::move((*header)[i]));
    }

    auto logger = gelex::logging::get();
    logger->info(
        "Loaded {} samples with {} qcovars", map->size(), covar_names.size());

    if (covar_names.size() <= 3)
    {
        logger->info("qcovar names: {}", fmt::join(covar_names, ", "));
    }
    else
    {
        logger->info(
            "qcovar names: {}, {}, {}, ... ({} total)",
            covar_names[0],
            covar_names[1],
            covar_names[2],
            covar_names.size());
    }

    return QcovarLoader(std::move(covar_names), std::move(*map));
}

void QcovarLoader::intersect(std::unordered_set<std::string>& id_set) const
{
    std::erase_if(
        id_set,
        [this](const std::string& current_id)
        { return !covariate_data_.contains(current_id); });
}

Eigen::MatrixXd QcovarLoader::load(
    const std::unordered_map<std::string, Eigen::Index>& id_map) const
{
    Eigen::MatrixXd data(id_map.size(), covariate_names_.size());
    data.setZero();

    for (const auto& [id, values] : covariate_data_)
    {
        if (auto it = id_map.find(id); it != id_map.end())
        {
            for (Eigen::Index i = 0; i < values.size(); ++i)
            {
                data(it->second, i) = values[i];
            }
        }
    }
    return data;
}

auto QcovarLoader::get_header(std::ifstream& file)
    -> std::expected<std::vector<std::string>, Error>
{
    std::string line;
    std::getline(file, line);
    const auto header_view = parse_header(line);

    if (!header_view)
    {
        return std::unexpected(header_view.error());
    }

    if (header_view->size() < 3)
    {
        return std::unexpected(
            Error{
                ErrorCode::InvalidRange,
                std::format(
                    "Qcovar file must have at least 3 columns, got {}",
                    header_view->size())});
    }
    return std::ranges::to<std::vector<std::string>>(*header_view);
}

auto QcovarLoader::read(
    std::ifstream& file,
    size_t expected_columns,
    bool iid_only) -> std::
    expected<std::unordered_map<std::string, std::vector<double>>, Error>
{
    std::unordered_map<std::string, std::vector<double>> covariate_data;

    std::string line;
    int n_line{};
    while (std::getline(file, line))
    {
        if (line.empty())
        {
            continue;
        }
        if (count_num_columns(line) != expected_columns)
        {
            return std::unexpected(
                Error{
                    ErrorCode::InconsistColumnCount,
                    std::format(
                        "Inconsistent number of columns (line {})",
                        n_line + 2)});
        }

        if (auto id_str = parse_id(line, iid_only); id_str)
        {
            auto values = parse_all_doubles(line, 2);  // skip FID and IID
            if (!values)
            {
                return std::unexpected(enrich_with_line_info(
                    std::move(values.error()), n_line + 2));
            }

            covariate_data.emplace(std::move(*id_str), std::move(*values));
        }
        else
        {
            return std::unexpected(
                enrich_with_line_info(std::move(id_str.error()), n_line + 2));
        }
        n_line++;
    }
    return covariate_data;
}

auto CovarLoader::create(std::string_view path, bool iid_only)
    -> std::expected<CovarLoader, Error>
{
    auto file = detail::openfile<std::ifstream>(path);
    if (!file)
    {
        return std::unexpected(file.error());
    }
    auto header = get_header(*file);
    if (!header)
    {
        return std::unexpected(
            enrich_with_file_info(std::move(header.error()), path));
    }
    size_t cols = header->size();
    auto map = read(*file, cols, iid_only);
    if (!map)
    {
        return std::unexpected(
            enrich_with_file_info(std::move(map.error()), path));
    }

    std::vector<std::string> covar_names;
    covar_names.reserve(cols - 2);
    for (size_t i = 2; i < cols; ++i)
    {
        covar_names.push_back(std::move((*header)[i]));
    }

    auto encode_maps = build_encode_maps(*map, covar_names);

    auto logger = gelex::logging::get();

    logger->info(
        "Loaded {} samples with {} covars", map->size(), covar_names.size());

    std::vector<std::string> covar_names_with_levels;
    for (size_t i = 0; i < covar_names.size(); ++i)
    {
        size_t num_levels = encode_maps[i].size();
        if (num_levels > 0)
        {
            covar_names_with_levels.push_back(
                std::format("{} [{}]", covar_names[i], num_levels));
        }
        else
        {
            covar_names_with_levels.push_back(covar_names[i]);
        }
    }

    if (covar_names_with_levels.size() <= 3)
    {
        logger->info(
            "covar names: {}", fmt::join(covar_names_with_levels, ", "));
    }
    else
    {
        logger->info(
            "covar names: {}, {}, {}, ... ({} total)",
            covar_names_with_levels[0],
            covar_names_with_levels[1],
            covar_names_with_levels[2],
            covar_names_with_levels.size());
    }

    return CovarLoader(
        std::move(covar_names), std::move(*map), std::move(encode_maps));
}

void CovarLoader::intersect(std::unordered_set<std::string>& id_set)
{
    size_t initial_count = id_set.size();

    std::erase_if(
        covariate_data_,
        [&id_set](const auto& pair) { return !id_set.contains(pair.first); });

    std::erase_if(
        id_set,
        [this](const std::string& current_id)
        { return !covariate_data_.contains(current_id); });

    rebuild_encode_maps();
}

Eigen::MatrixXd CovarLoader::load(
    const std::unordered_map<std::string, Eigen::Index>& id_map) const
{
    size_t total_dummy_vars = 0;
    std::vector<Eigen::Index> col_offsets;
    for (const auto& encode_map : encode_maps_)
    {
        col_offsets.push_back(static_cast<Eigen::Index>(total_dummy_vars));
        if (!encode_map.empty())
        {
            total_dummy_vars += encode_map.begin()->second.size();
        }
    }

    Eigen::MatrixXd data(id_map.size(), total_dummy_vars);
    data.setZero();

    for (const auto& [id, values] : covariate_data_)
    {
        const auto id_it = id_map.find(id);
        if (id_it == id_map.end())
        {
            continue;
        }
        const Eigen::Index row_idx = id_it->second;

        const size_t num_covars_to_process
            = std::min(values.size(), encode_maps_.size());
        for (size_t i = 0; i < num_covars_to_process; ++i)
        {
            if (const auto map_it = encode_maps_[i].find(values[i]);
                map_it != encode_maps_[i].end())
            {
                const auto& dummy_encoding = map_it->second;
                const Eigen::Index offset = col_offsets[i];
                const auto num_vars
                    = static_cast<Eigen::Index>(dummy_encoding.size());

                data.row(row_idx).segment(offset, num_vars)
                    = Eigen::Map<const Eigen::VectorXi>(
                          dummy_encoding.data(), num_vars)
                          .template cast<double>();
            }
        }
    }
    return data;
}

auto CovarLoader::get_header(std::ifstream& file)
    -> std::expected<std::vector<std::string>, Error>
{
    std::string line;
    std::getline(file, line);
    const auto header_view = parse_header(line);

    if (!header_view)
    {
        return std::unexpected(header_view.error());
    }

    if (header_view->size() < 3)
    {
        return std::unexpected(
            Error{
                ErrorCode::InvalidRange,
                std::format(
                    "Covar file must have at least 3 columns, got {}",
                    header_view->size())});
    }
    return std::ranges::to<std::vector<std::string>>(*header_view);
}

auto CovarLoader::read(
    std::ifstream& file,
    size_t expected_columns,
    bool iid_only) -> std::
    expected<std::unordered_map<std::string, std::vector<std::string>>, Error>
{
    std::unordered_map<std::string, std::vector<std::string>> covariate_data;

    std::string line;
    int n_line{};
    while (std::getline(file, line))
    {
        if (line.empty())
        {
            continue;
        }
        if (count_num_columns(line) != expected_columns)
        {
            return std::unexpected(
                Error{
                    ErrorCode::InconsistColumnCount,
                    std::format(
                        "Inconsistent number of columns at line {}",
                        n_line + 2)});
        }

        if (auto id_str = parse_id(line, iid_only); id_str)
        {
            auto values = parse_string(line, 2);  // skip FID and IID
            if (values.size() != expected_columns - 2)
            {
                return std::unexpected(
                    Error{
                        ErrorCode::InconsistColumnCount,
                        std::format(
                            "Inconsistent number of columns at line {}",
                            n_line + 2)});
            }

            std::vector<std::string> string_values;
            string_values.reserve(values.size());
            for (const auto& value : values)
            {
                string_values.emplace_back(value);
            }

            covariate_data.emplace(
                std::move(*id_str), std::move(string_values));
        }
        else
        {
            return std::unexpected(
                enrich_with_line_info(std::move(id_str.error()), n_line + 2));
        }
        n_line++;
    }
    return covariate_data;
}

auto CovarLoader::build_encode_maps(
    const std::unordered_map<std::string, std::vector<std::string>>&
        covariate_data,
    std::span<const std::string> covariate_names)
    -> std::vector<std::unordered_map<std::string, std::vector<int>>>
{
    if (covariate_names.empty())
    {
        return {};
    }

    auto unique_levels_per_covariate
        = collect_unique_levels(covariate_data, covariate_names.size());

    std::vector<std::unordered_map<std::string, std::vector<int>>> result;
    result.reserve(covariate_names.size());

    for (const auto& unique_levels : unique_levels_per_covariate)
    {
        result.push_back(create_encoding_for_one_covariate(unique_levels));
    }

    return result;
}

auto CovarLoader::collect_unique_levels(
    const std::unordered_map<std::string, std::vector<std::string>>& data,
    size_t num_covariates) -> std::vector<std::unordered_set<std::string_view>>
{
    std::vector<std::unordered_set<std::string_view>> all_levels(
        num_covariates);
    for (const auto& [id, values] : data)
    {
        for (size_t i = 0; i < num_covariates; ++i)
        {
            if (!values[i].empty())
            {
                all_levels[i].insert(values[i]);
            }
        }
    }
    return all_levels;
}

auto CovarLoader::create_encoding_for_one_covariate(
    const std::unordered_set<std::string_view>& unique_levels)
    -> std::unordered_map<std::string, std::vector<int>>
{
    if (unique_levels.size() < 2)
    {
        return {};
    }

    std::vector<std::string_view> sorted_levels(
        unique_levels.begin(), unique_levels.end());
    std::ranges::sort(sorted_levels);

    const size_t num_dummy_vars = sorted_levels.size() - 1;
    std::unordered_map<std::string, std::vector<int>> encode_map;
    encode_map.reserve(sorted_levels.size());

    encode_map.emplace(sorted_levels[0], std::vector<int>(num_dummy_vars, 0));

    for (size_t i = 1; i < sorted_levels.size(); ++i)
    {
        std::vector<int> encoding(num_dummy_vars, 0);
        encoding[i - 1] = 1;
        encode_map.emplace(sorted_levels[i], std::move(encoding));
    }
    return encode_map;
}

void CovarLoader::rebuild_encode_maps()
{
    encode_maps_ = build_encode_maps(covariate_data_, covariate_names_);
}

auto BimLoader::create(std::string_view path) -> std::expected<BimLoader, Error>
{
    auto file = detail::openfile<std::ifstream>(path);
    if (!file)
    {
        return std::unexpected(file.error());
    }
    auto snp_ids = read(*file);
    if (!snp_ids)
    {
        return std::unexpected(
            enrich_with_file_info(std::move(snp_ids.error()), path));
    }

    auto logger = gelex::logging::get();
    logger->info("Loaded {} SNP IDs from bim file", snp_ids->size());

    return BimLoader(std::move(*snp_ids));
}

auto BimLoader::read(std::ifstream& file)
    -> std::expected<std::vector<std::string>, Error>
{
    constexpr static size_t bim_n_cols = 6;
    std::vector<std::string> snp_ids;
    std::string line;
    int n_line{};
    while (std::getline(file, line))
    {
        if (line.empty())
        {
            continue;
        }

        if (count_num_columns(line) != bim_n_cols)
        {
            return std::unexpected(
                Error{
                    ErrorCode::InconsistColumnCount,
                    std::format(
                        "Bim file must have at least 2 columns at line {}",
                        n_line + 1)});
        }

        auto tokens = parse_string(line, 0, "\t");
        if (tokens.size() < 2)
        {
            return std::unexpected(
                Error{
                    ErrorCode::InvalidFile,
                    std::format(
                        "Failed to parse SNP ID from line {}", n_line + 1)});
        }

        snp_ids.emplace_back(tokens[1]);
        n_line++;
    }
    return snp_ids;
}

auto FamLoader::create(std::string_view path, bool iid_only)
    -> std::expected<FamLoader, Error>
{
    auto file = detail::openfile<std::ifstream>(path);
    if (!file)
    {
        return std::unexpected(file.error());
    }
    auto individual_ids = read(*file, iid_only);
    if (!individual_ids)
    {
        return std::unexpected(
            enrich_with_file_info(std::move(individual_ids.error()), path));
    }

    auto logger = gelex::logging::get();
    logger->info(
        "Loaded {} samples from fam file[{}]", individual_ids->size(), path);

    return FamLoader(std::move(*individual_ids));
}

void FamLoader::intersect(std::unordered_set<std::string>& id_set) const
{
    std::erase_if(
        id_set,
        [this](const std::string& current_id)
        { return !sample_ids_.contains(current_id); });
}

auto FamLoader::read(std::ifstream& file, bool iid_only)
    -> std::expected<std::unordered_set<std::string>, Error>
{
    constexpr static size_t fam_n_cols = 6;
    std::unordered_set<std::string> individual_ids;
    std::string line;
    int n_line{};
    while (std::getline(file, line))
    {
        if (line.empty())
        {
            continue;
        }

        if (count_num_columns(line, " ") != fam_n_cols)
        {
            return std::unexpected(
                Error{
                    ErrorCode::InconsistColumnCount,
                    std::format(
                        ".fam file must have at 6 columns (line {})",
                        n_line + 1)});
        }

        if (auto id_str = parse_id(line, iid_only, " "); id_str)
        {
            individual_ids.insert(std::move(*id_str));
        }
        else
        {
            return std::unexpected(
                enrich_with_line_info(std::move(id_str.error()), n_line + 1));
        }
        n_line++;
    }
    return individual_ids;
}

}  // namespace gelex::detail
