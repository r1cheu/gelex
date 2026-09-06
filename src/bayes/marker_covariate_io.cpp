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

#include <algorithm>
#include <array>
#include <filesystem>
#include <ranges>
#include <string>
#include <string_view>

#include "gelex/data/dataframe/column.h"
#include "gelex/data/dataframe/dataframe.h"
#include "gelex/data/dataframe/reader.h"
#include "gelex/exception.h"

namespace gelex::bayes
{
auto read_marker_annotation(const std::filesystem::path& path)
    -> DataFrame<std::string>
{
    using enum ColumnType;
    constexpr std::array schema{String, Int, String, String, Double};
    constexpr std::array<std::string_view, 4> metadata_names{
        "CHR", "BP", "A1", "A2"};

    ReadOptions options;
    options.delimiter = '\t';
    options.index_cols = {1};
    auto frame = read_dataframe<std::string>(path, options, schema);
    if (!std::ranges::equal(
            frame.names() | std::views::take(metadata_names.size()),
            metadata_names))
    {
        throw GelexException(
            "marker annotation header must be "
            "CHR\tSNP\tBP\tA1\tA2\t<annotation>");
    }
    return frame;
}
}  // namespace gelex::bayes
