#include "gelex/data/io.h"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstddef>
#include <expected>
#include <format>
#include <fstream>
#include <iterator>
#include <ranges>
#include <stdexcept>
#include <string>
#include <string_view>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

#include <fmt/format.h>
#include <fmt/ranges.h>
#include <Eigen/Core>

#include "../src/data/parser.h"

namespace gelex
{
namespace detail
{

void validate_path_or_throw(const std::filesystem::path& path)
{
    if (!path.empty() && !std::filesystem::exists(path))
    {
        throw std::runtime_error(
            std::format("Required file does not exist: {}", path.string()));
    }
}
}  // namespace detail

using Eigen::Index;

using Eigen::MatrixXd;
using Eigen::VectorXd;

DataReader DataReader::Create(
    const std::filesystem::path& pheno_path,
    const std::filesystem::path& fam_path,
    const std::filesystem::path& qcovar_path,
    const std::filesystem::path& covar_path,
    size_t pheno_col,
    bool iid_only)
{
    // check must exist
    detail::validate_path_or_throw(pheno_path);
    detail::validate_path_or_throw(fam_path);
    bool has_qcovar = !qcovar_path.empty();
    bool has_covar = !covar_path.empty();

    // parse phenotype and intersect
    id_set phenotype_ids = parse_phenotype(pheno_path, pheno_col, iid_only);
    std::unordered_set<std::string> interseted_ids = phenotype_ids;
    id_set fam_ids = parse_fam(fam_path, iid_only);
    intersect_in_place(interseted_ids, fam_ids);

    // if qcovar is provided, intersect
    if (has_qcovar)
    {
        detail::validate_path_or_throw(qcovar_path);
        id_set qcovar_ids = parse_covar(qcovar_path, iid_only);
        intersect_in_place(interseted_ids, qcovar_ids);
    }

    // if covar is provided, intersect
    if (has_covar)
    {
        detail::validate_path_or_throw(covar_path);
        id_set covar_ids = parse_covar(covar_path, iid_only);
        intersect_in_place(interseted_ids, covar_ids);
    }

    // move the set to a sorted vector
    std::vector<std::string> final_ids(
        std::make_move_iterator(interseted_ids.begin()),
        std::make_move_iterator(interseted_ids.end()));
    std::ranges::sort(final_ids);

    // create id to index map
    std::unordered_map<std::string, Eigen::Index> id_map;
    for (Index i = 0; i < final_ids.size(); ++i)
    {
        id_map[final_ids[i]] = i;
    }
    auto n_individuals = static_cast<Index>(id_map.size());

    VectorXd phenotype
        = load_phenotype(pheno_path, id_map, iid_only, pheno_col);

    MatrixXd qcovar;
    if (has_qcovar)
    {
        qcovar = std::move(load_qcovar(qcovar_path, id_map, iid_only));
    }

    MatrixXd covar;
    if (has_covar)
    {
        covar = std::move(load_covar(covar_path, id_map, iid_only));
    }

    Index n_fixed = 1 + qcovar.cols() + covar.cols();

    MatrixXd fixed(n_individuals, n_fixed);
    fixed.col(0) = VectorXd::Ones(n_individuals);
    fixed.middleCols(1, qcovar.cols()) = qcovar;
    fixed.rightCols(covar.cols()) = covar;

    return DataReader(
        std::move(phenotype),
        std::move(fixed),
        std::move(final_ids),
        std::move(id_map));
}

VectorXd DataReader::load_phenotype(
    const std::filesystem::path& pheno_path,
    const std::unordered_map<std::string, Eigen::Index>& id_map,
    bool iid_only,
    size_t pheno_col)
{
    // resize phenotype vector
    VectorXd phenotype(id_map.size());

    auto file = detail::open_or_throw<std::ifstream>(pheno_path);
    std::string line;
    std::getline(file, line);

    while (std::getline(file, line))
    {
        if (line.empty())
        {
            continue;
        }

        // we don't check id parsing error here, because the invalid id should
        // be handle in parsing step
        if (auto id_str = parse_id(line, iid_only); id_str)
        {
            auto it = id_map.find(*id_str);
            if (it != id_map.end())
            {
                phenotype(it->second)
                    = detail::parse_nth_double(line, pheno_col - 1).value();
            }
        }
    }
    return phenotype;
}

MatrixXd DataReader::load_qcovar(
    const std::filesystem::path& qcovar_path,
    const std::unordered_map<std::string, Eigen::Index>& id_map,
    bool iid_only)
{
    auto file = detail::open_or_throw<std::ifstream>(qcovar_path);

    std::string line;
    std::getline(file, line);
    auto header = parse_header(line, qcovar_path);

    auto rows = static_cast<Index>(id_map.size());
    auto cols = static_cast<Index>(header.size() - 2);

    MatrixXd qcovar(rows, cols);

    while (std::getline(file, line))
    {
        if (line.empty())
        {
            continue;
        }
        if (auto id_str = parse_id(line, iid_only); id_str)
        {
            auto it = id_map.find(*id_str);
            if (it != id_map.end())
            {
                auto qcovar_i = detail::parse_all_doubles(line, 2).value();
                qcovar.row(it->second) = Eigen::Map<Eigen::RowVectorXd>(
                    qcovar_i.data(), static_cast<Index>(qcovar_i.size()));
            }
        }
    }
    return qcovar;
}

MatrixXd DataReader::load_covar(
    const std::filesystem::path& covar_path,
    const std::unordered_map<std::string, Eigen::Index>& id_map,
    bool iid_only)
{
    covariate_encoding_map encoding_map
        = encode_covar(covar_path, id_map, iid_only);
    auto file = detail::open_or_throw<std::ifstream>(covar_path);

    std::string line;
    std::getline(file, line);

    auto rows = static_cast<Index>(id_map.size());
    Index cols = 0;
    for (const auto& covar : encoding_map)
    {
        cols += static_cast<Index>(covar.second.size() - 1);
    }

    MatrixXd covar(rows, cols);

    while (std::getline(file, line))
    {
        if (line.empty())
        {
            continue;
        }

        if (auto id_str = parse_id(line, iid_only); id_str)
        {
            auto it = id_map.find(*id_str);
            if (it != id_map.end())
            {
                std::vector<std::string_view> qcovar_i
                    = detail::parse_string(line, 2);
                std::vector<double> row_i;
                row_i.reserve(cols);
                for (auto&& [level, level_encode] :
                     std::views::zip(qcovar_i, encoding_map))
                {
                    auto level_it
                        = level_encode.second.find(std::string(level));
                    if (level_it != level_encode.second.end())
                    {
                        row_i.insert(
                            row_i.end(),
                            level_it->second.begin(),
                            level_it->second.end());
                    }
                }
                covar.row(it->second) = Eigen::Map<Eigen::RowVectorXd>(
                    row_i.data(), static_cast<Index>(row_i.size()));
            }
        }
    }
    return covar;
}

DataReader::covariate_encoding_map DataReader::encode_covar(
    const std::filesystem::path& path,
    const std::unordered_map<std::string, Eigen::Index>& id_map,
    bool iid_only)
{
    auto file = detail::open_or_throw<std::ifstream>(path);

    std::string line;
    std::getline(file, line);
    const auto covar_name = detail::parse_string(line, 2)
                            | std::ranges::to<std::vector<std::string>>();
    const size_t n_covar = covar_name.size();

    std::vector<std::unordered_set<std::string>> unique_values(n_covar);

    size_t n_line = 2;
    while (std::getline(file, line))
    {
        if (line.empty())
        {
            continue;
        }
        auto id_str = parse_id(line, iid_only);
        if (id_str && id_map.contains(*id_str))
        {
            const auto values = detail::parse_string(line, 2);
            if (values.size() != n_covar)
            {
                throw std::runtime_error(
                    std::format(
                        "Number of columns in covar file [{}] line {} is not "
                        "consist.",
                        path.string(),
                        n_line++));
            }
            auto value_it = values.begin();
            for (size_t i = 0; i < n_covar && value_it != values.end();
                 ++i, ++value_it)
            {
                unique_values[i].emplace(*value_it);
            }
        }
    }

    covariate_encoding_map final_encoding_map;
    final_encoding_map.reserve(n_covar);
    for (auto&& [name, values_set] : std::views::zip(covar_name, unique_values))
    {
        if (values_set.size() < 2)
        {
            continue;
        }

        std::vector<std::string> sorted_levels(
            std::make_move_iterator(values_set.begin()),
            std::make_move_iterator(values_set.end()));
        std::ranges::sort(sorted_levels);

        const size_t encoding_dim = sorted_levels.size() - 1;
        level_encode_map level_map;
        level_map.reserve(sorted_levels.size());

        level_map.emplace(sorted_levels[0], encode_vector(encoding_dim, 0.0));

        for (size_t i = 1; i < sorted_levels.size(); ++i)
        {
            encode_vector one_hot(encoding_dim, 0.0);
            one_hot[i - 1] = 1.0;
            level_map.emplace(std::move(sorted_levels[i]), std::move(one_hot));
        }

        final_encoding_map.emplace(std::string(name), std::move(level_map));
    }

    return final_encoding_map;
}

DataReader::id_set DataReader::parse_phenotype(
    const std::filesystem::path& path,
    size_t pheno_col,
    bool iid_only)
{
    if (pheno_col < 3)
    {
        throw std::invalid_argument(
            std::format(
                "Phenotype column must be 3 or greater, but got {}",
                pheno_col));
    }

    auto file = detail::open_or_throw<std::ifstream>(path);
    id_set result;
    std::string line;
    std::getline(file, line);

    auto header = parse_header(line, path);

    if (header.size() < pheno_col)
    {
        throw std::invalid_argument(
            std::format(
                "Phenotyhpe file [{}] has only {} columns, but column {} was "
                "requested",
                path.string(),
                header.size(),
                pheno_col));
    }

    size_t n_line = 2;
    while (std::getline(file, line))
    {
        if (line.empty())
        {
            continue;
        }
        if (auto id_str = parse_id(line, iid_only); id_str)
        {
            if (auto val = detail::parse_nth_double(line, pheno_col - 1); val)
            {
                result.emplace(std::move(*id_str));
            }
        }
        else
        {
            throw std::runtime_error(
                std::format(
                    "Failed to parse id of line {:d} in Phenotype file [{}]",
                    n_line++,
                    path.string()));
        }
    }
    return result;
};

DataReader::id_set DataReader::parse_covar(
    const std::filesystem::path& path,
    bool iid_only)
{
    auto file = detail::open_or_throw<std::ifstream>(path);
    id_set result;
    std::string line;
    std::getline(file, line);
    constexpr size_t id_cols = 2;

    auto header = parse_header(line, path);

    size_t n_line = 2;  // skip header
    while (std::getline(file, line))
    {
        if (line.empty())
        {
            continue;
        }
        if (auto id_str = parse_id(line, iid_only); id_str)
        {
            result.emplace(std::move(*id_str));
        }
        else
        {
            throw std::runtime_error(
                std::format(
                    "Failed to parse id of line {:d} in qcovar file [{}]",
                    n_line++,
                    path.string()));
        }
    }
    return result;
};

DataReader::id_set DataReader::parse_fam(
    const std::filesystem::path& path,
    bool iid_only)
{
    auto file = detail::open_or_throw<std::ifstream>(path);
    id_set result;
    std::string line;

    size_t n_line = 2;
    while (std::getline(file, line))
    {
        if (line.empty())
        {
            continue;
        }
        if (auto id_str = parse_id(line, iid_only, " "); id_str)
        {
            result.emplace(std::move(*id_str));
        }
        else
        {
            throw std::runtime_error(
                std::format(
                    "Failed to parse id of line {:d} in fam file [{}]",
                    n_line++,
                    path.string()));
        }
    }
    return result;
};

void DataReader::intersect_in_place(
    std::unordered_set<std::string>& main_set,
    const std::unordered_set<std::string>& other_set)
{
    for (auto it = main_set.begin(); it != main_set.end();)
    {
        if (!other_set.contains(*it))
        {
            it = main_set.erase(it);
        }
        else
        {
            ++it;
        }
    }
}
}  // namespace gelex
