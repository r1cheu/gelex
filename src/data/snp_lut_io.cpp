// Copyright 2026 RuLei Chen
// SPDX-License-Identifier: Apache-2.0

#include "gelex/data/snp_lut_io.h"

#include <Eigen/Core>
#include <cstddef>
#include <cstdint>
#include <filesystem>
#include <fmt/format.h>

#include "gelex/exception.h"
#include "gelex/genetic_mode.h"
#include "gelex/io/binary_reader.h"
#include "gelex/io/dense_writer.h"

namespace gelex
{

auto load_snp_luts(const std::filesystem::path& path) -> ModeMap<SnpLutMatrix>
{
    BinaryReader reader(path.string());
    ModeMap<SnpLutMatrix> luts;
    for (const auto mode : all_genetic_modes)
    {
        if (reader.contains(fmt::format("{}/lut", mode)))
        {
            auto lut_map = reader.to_map<double>(fmt::format("{}/lut", mode));
            if (lut_map.rows() != 4)
            {
                throw GelexException(
                    fmt::format(
                        "load_snp_luts: mode {} requires 4 LUT rows, got {}",
                        mode,
                        lut_map.rows()));
            }
            luts.emplace(mode, lut_map.array());
        }
    }
    return luts;
}

auto write_snp_luts(
    const std::filesystem::path& path,
    const ModeMap<SnpLutMatrix>& luts) -> void
{
    auto writer = open_dense_writer(path.string());
    for (const auto& [mode, lut] : luts)
    {
        writer.reserve<double>(
            fmt::format("{}/lut", mode),
            BinaryShape{
                static_cast<std::uint64_t>(lut.rows()),
                static_cast<std::uint64_t>(lut.cols())})
            << lut.reshaped();
    }
    writer.close();
}

}  // namespace gelex
