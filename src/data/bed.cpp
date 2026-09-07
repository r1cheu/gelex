// Copyright 2026 RuLei Chen
// SPDX-License-Identifier: Apache-2.0

#include "gelex/data/bed.h"

#include <filesystem>
#include <string>
#include <utility>

#include "gelex/data/detail/bed_source.h"
#include "gelex/data/reader.h"

namespace gelex
{

auto open_bed(const std::string& bfile_prefix) -> Bed
{
    auto source_index
        = read_fam(std::filesystem::path{bfile_prefix + ".fam"}).index();
    auto bim = read_bim(std::filesystem::path{bfile_prefix + ".bim"});
    auto bed_source = detail::open_bed_source(bfile_prefix);

    return Bed{std::move(bed_source), std::move(source_index), std::move(bim)};
}

}  // namespace gelex
