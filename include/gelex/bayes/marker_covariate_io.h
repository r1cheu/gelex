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

#ifndef GELEX_BAYES_MARKER_COVARIATE_IO_H
#define GELEX_BAYES_MARKER_COVARIATE_IO_H

#include <filesystem>
#include <string>

#include "gelex/data/dataframe/key_type.h"

namespace gelex
{
template <KeyType Key>
class DataFrame;
}  // namespace gelex

namespace gelex::bayes
{
// Reads a marker annotation table with header CHR\tSNP\tBP\tA1\tA2\t<name>
// (one numeric annotation column), indexed by SNP. Alignment against marker
// metadata is make_marker_covariate.
auto read_marker_annotation(const std::filesystem::path& path)
    -> DataFrame<std::string>;
}  // namespace gelex::bayes

#endif  // GELEX_BAYES_MARKER_COVARIATE_IO_H
