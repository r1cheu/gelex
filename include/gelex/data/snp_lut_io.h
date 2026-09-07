// Copyright 2026 RuLei Chen
// SPDX-License-Identifier: Apache-2.0

#ifndef GELEX_DATA_SNP_LUT_IO_H_
#define GELEX_DATA_SNP_LUT_IO_H_

#include <filesystem>

#include "gelex/data/snp_lut.h"
#include "gelex/genetic_mode.h"

namespace gelex
{

[[nodiscard]] auto load_snp_luts(const std::filesystem::path& path)
    -> ModeMap<SnpLutMatrix>;

auto write_snp_luts(
    const std::filesystem::path& path,
    const ModeMap<SnpLutMatrix>& luts) -> void;

}  // namespace gelex

#endif  // GELEX_DATA_SNP_LUT_IO_H_
