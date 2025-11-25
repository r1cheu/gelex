#include "covariate_loader.h"

#include <filesystem>
#include <fstream>
#include <map>
#include <ranges>
#include <string>
#include <unordered_map>
#include <vector>

#include <fmt/ranges.h>
#include <Eigen/Core>

#include "../src/data/loader.h"
#include "../src/data/parser.h"
#include "gelex/exception.h"
#include "gelex/logger.h"

namespace gelex::detail
{

CovarPredictLoader::CovarPredictLoader(
    const std::filesystem::path& path,
    bool iid_only)
{
    auto file = detail::open_file<std::ifstream>(path, std::ios_base::in);

    auto header = get_header(file);
    size_t cols = header.size();
    auto map = read(file, cols, iid_only);

    std::vector<std::string> covar_names;
    covar_names.reserve(cols - 2);
    for (size_t i = 2; i < cols; ++i)
    {
        covar_names.push_back(std::move(header[i]));
    }

    names_ = std::move(covar_names);
    data_ = std::move(map);

    auto logger = gelex::logging::get();
    logger->info(
        "Loaded {} samples with {} covars.", data_.size(), names_.size());
}

auto CovarPredictLoader::load(
    const std::unordered_map<std::string, Eigen::Index>& id_map) const
    -> std::map<std::string, std::vector<std::string>>
{
    const auto n_samples = id_map.size();

    std::map<std::string, std::vector<std::string>> formatted_data;

    std::vector<std::vector<std::string>*> col_ptrs;
    col_ptrs.reserve(names_.size());

    for (const auto& covar_name : names_)
    {
        formatted_data[covar_name] = std::vector<std::string>(n_samples);
        col_ptrs.push_back(&formatted_data[covar_name]);
    }

    for (const auto& [id, row_idx] : id_map)
    {
        auto data_it = data_.find(id);
        if (data_it != data_.end())
        {
            const auto& values = data_it->second;

            if (values.size() != col_ptrs.size())
            {
                continue;
            }

            for (size_t i = 0; i < col_ptrs.size(); ++i)
            {
                (*col_ptrs[i])[row_idx] = values[i];
            }
        }
    }
    return formatted_data;
}

auto CovarPredictLoader::get_header(std::ifstream& file)
    -> std::vector<std::string>
{
    std::string line;
    std::getline(file, line);
    const auto header_view = parse_header(line);

    if (header_view.size() < 3)
    {
        throw InvalidRangeException(
            std::format(
                "Covar file must have at least 3 columns, got {}",
                header_view.size()));
    }
    return std::ranges::to<std::vector<std::string>>(header_view);
}

auto CovarPredictLoader::read(
    std::ifstream& file,
    size_t expected_columns,
    bool iid_only) -> std::unordered_map<std::string, std::vector<std::string>>
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
            throw InconsistentColumnCountException(
                std::format(
                    "Inconsistent number of columns at line {}", n_line + 2));
        }

        auto id_str = parse_id(line, iid_only);
        auto values = parse_string(line, 2);  // skip FID and IID
        if (values.size() != expected_columns - 2)
        {
            throw InconsistentColumnCountException(
                std::format(
                    "Inconsistent number of columns at line {}", n_line + 2));
        }

        std::vector<std::string> string_values;
        string_values.reserve(values.size());
        for (const auto& value : values)
        {
            string_values.emplace_back(value);
        }

        covariate_data.emplace(std::move(id_str), std::move(string_values));
        n_line++;
    }
    return covariate_data;
}
}  // namespace gelex::detail
