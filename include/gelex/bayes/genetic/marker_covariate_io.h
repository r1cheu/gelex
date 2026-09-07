// Copyright 2026 RuLei Chen
// SPDX-License-Identifier: Apache-2.0

#ifndef GELEX_BAYES_GENETIC_MARKER_COVARIATE_IO_H
#define GELEX_BAYES_GENETIC_MARKER_COVARIATE_IO_H

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

#endif  // GELEX_BAYES_GENETIC_MARKER_COVARIATE_IO_H
