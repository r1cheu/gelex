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

#include "fam_loader.h"

#include <filesystem>
#include <fstream>
#include <ranges>
#include <vector>

#include <Eigen/Core>

#include <fmt/ranges.h>
#include "../parser.h"
#include "gelex/exception.h"

namespace gelex::detail
{

FamLoader::FamLoader(const std::filesystem::path& path, bool iid_only)
{
    try
    {
        set_ids(path, iid_only);
        set_index_map();
    }
    catch (const GelexException& e)
    {
        throw FileFormatException(
            std::format("{}:{}", path.string(), e.what()));
    }
}

void FamLoader::set_ids(const std::filesystem::path& path, bool iid_only)
{
    ids_.clear();
    auto file = detail::open_file<std::ifstream>(path, std::ios::in);
    char delimiter = detect_file_delimiter(file);

    ids_.reserve(1024);  // Start small
    std::string line;
    int n_line = 0;

    while (std::getline(file, line))
    {
        n_line++;
        if (line.empty())
        {
            continue;
        }
        try
        {
            ids_.emplace_back(parse_id(line, iid_only, delimiter));
        }
        catch (const GelexException& e)
        {
            throw DataParseException(std::format("{}: {}", n_line, e.what()));
        }
    }
}

void FamLoader::set_index_map()
{
    data_.reserve(ids_.size());
    for (auto&& [idx, id] : ids_ | std::views::enumerate)
    {
        data_[id] = static_cast<Eigen::Index>(idx);
    }
}

}  // namespace gelex::detail
