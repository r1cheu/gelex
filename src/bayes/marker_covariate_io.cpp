// Copyright 2026 RuLei Chen
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//     http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.

#include "gelex/bayes/marker_covariate_io.h"

#include <array>
#include <cstddef>
#include <filesystem>
#include <fstream>
#include <ranges>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

#include "gelex/data/dataframe/column.h"
#include "gelex/data/dataframe/dataframe.h"
#include "gelex/data/dataframe/reader.h"
#include "gelex/exception.h"
#include "gelex/io/detail/parser.h"

namespace gelex::bayes
{
auto read_marker_covariate(
    const std::filesystem::path& path,
    const DataFrame<std::string>& marker_metadata) -> MarkerCovariate
{
    auto stream = detail::open_file<std::ifstream>(path, std::ios::in);
    std::string header;
    if (!std::getline(stream, header))
    {
        throw GelexException("marker covariate file is empty");
    }
    if (!header.empty() && header.back() == '\r')
    {
        header.pop_back();
    }

    std::vector<std::string> names;
    std::size_t begin = 0;
    while (true)
    {
        const auto end = header.find('\t', begin);
        names.emplace_back(header.substr(begin, end - begin));
        if (end == std::string::npos)
        {
            break;
        }
        begin = end + 1;
    }

    constexpr std::array<std::string_view, 5> required_names{
        "CHR", "SNP", "BP", "A1", "A2"};
    if (names.size() <= required_names.size())
    {
        throw GelexException(
            "marker covariate file requires at least one annotation column");
    }
    for (const auto& [name, required_name] :
         std::views::zip(names, required_names))
    {
        if (name != required_name)
        {
            throw GelexException(
                "marker covariate header must begin with "
                "CHR\tSNP\tBP\tA1\tA2");
        }
    }

    std::vector<ColumnType> schema{
        ColumnType::String,
        ColumnType::Int,
        ColumnType::String,
        ColumnType::String};
    schema.resize(names.size() - 1, ColumnType::Double);

    ReadOptions options;
    options.delimiter = '\t';
    options.index_cols = {1};
    auto frame = read_dataframe<std::string>(path, options, schema);
    return make_marker_covariate(std::move(frame), marker_metadata);
}
}  // namespace gelex::bayes
