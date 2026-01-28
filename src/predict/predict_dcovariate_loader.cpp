/*
 * Copyright 2026 RuLei Chen
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 * http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

#include "predict_dcovariate_loader.h"

#include <filesystem>
#include <fstream>
#include <map>
#include <ranges>
#include <string>
#include <unordered_map>
#include <vector>

#include <fmt/ranges.h>
#include <Eigen/Core>

#include "../src/data/parser.h"
#include "gelex/exception.h"

namespace gelex::detail
{

DcovarPredictLoader::DcovarPredictLoader(
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
}

auto DcovarPredictLoader::load(
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

void DcovarPredictLoader::set_names(std::ifstream& file)
{
    std::string line;
    std::getline(file, line);
    const auto header_view = parse_header(line);

    if (header_view.size() < 3)
    {
        throw ColumnRangeException(
            std::format(
                "Covar file must have at least 3 columns, got {}",
                header_view.size()));
    }
    names_.clear();
    for (size_t i = 2; i < header_view.size(); ++i)
    {
        names_.emplace_back(header_view[i]);
    }
}

void DcovarPredictLoader::set_data(std::ifstream& file, bool iid_only)
{
    const size_t expected_columns = names_.size() + 2;
    data_.clear();
    data_.reserve(1024);

    std::string line;
    int n_line{};
    std::vector<std::string_view> value_buffer;
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
        parse_string(line, value_buffer, 2);  // skip FID and IID
        if (value_buffer.size() != expected_columns - 2)
        {
            throw InconsistentColumnCountException(
                std::format(
                    "Inconsistent number of columns at line {}", n_line + 2));
        }

        data_.emplace(
            std::move(id_str),
            std::ranges::to<std::vector<std::string>>(value_buffer));
        n_line++;
    }
}
}  // namespace gelex::detail
