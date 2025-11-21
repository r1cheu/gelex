#include "covariate_loader.h"

#include <expected>
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
#include "gelex/error.h"
#include "gelex/logger.h"

namespace gelex::detail
{

auto CovarPredictLoader::create(
    const std::filesystem::path& path,
    bool iid_only) -> std::expected<CovarPredictLoader, Error>
{
    auto file = detail::open_file<std::ifstream>(path, std::ios_base::in);
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
        "Loaded {} samples with {} covars.", map->size(), covar_names.size());

    return CovarPredictLoader(std::move(covar_names), std::move(*map));
}

auto CovarPredictLoader::load(
    const std::unordered_map<std::string, Eigen::Index>& id_map) const
    -> std::map<std::string, std::vector<std::string>>
{
    auto covariate_data = data_;
    std::erase_if(
        covariate_data,
        [&id_map](const auto& pair) { return !id_map.contains(pair.first); });

    return covariate_data;
}

auto CovarPredictLoader::get_header(std::ifstream& file)
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

auto CovarPredictLoader::read(
    std::ifstream& file,
    size_t expected_columns,
    bool iid_only)
    -> std::expected<std::map<std::string, std::vector<std::string>>, Error>
{
    std::map<std::string, std::vector<std::string>> covariate_data;

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
}  // namespace gelex::detail
