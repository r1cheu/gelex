// Copyright 2026 RuLei Chen
// SPDX-License-Identifier: Apache-2.0

#ifndef GELEX_DATA_SAMPLE_ID_IO_H_
#define GELEX_DATA_SAMPLE_ID_IO_H_

#include <span>
#include <string>

#include "gelex/data/dataframe/index.h"

namespace gelex
{

// Sample IDs travel as <prefix>.id: one "FID<TAB>IID" per line without a
// header, the layout PLINK's --keep accepts. ids are canonical FID<US>IID keys.
auto write_sample_ids(
    const std::string& prefix,
    std::span<const std::string> ids) -> void;

// Throws GelexException when <prefix>.id holds no samples.
[[nodiscard]] auto read_sample_ids(const std::string& prefix)
    -> DataFrameIndex<std::string>;

}  // namespace gelex

#endif  // GELEX_DATA_SAMPLE_ID_IO_H_
