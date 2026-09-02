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

#include "gelex/data/snp_lut_io.h"

#include <Eigen/Core>
#include <cstddef>
#include <cstdint>
#include <filesystem>
#include <fmt/format.h>
#include <span>

#include "gelex/exception.h"
#include "gelex/genetic_mode.h"
#include "gelex/io/binary_reader.h"
#include "gelex/io/binary_writer.h"

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
    BinaryWriter writer(path.string());
    for (const auto& [mode, lut] : luts)
    {
        writer
            .reserve<double>(
                fmt::format("{}/lut", mode),
                BinaryShape{
                    static_cast<std::uint64_t>(lut.rows()),
                    static_cast<std::uint64_t>(lut.cols())})
            .write(
                std::span<const double>{
                    lut.data(), static_cast<std::size_t>(lut.size())});
    }
}

}  // namespace gelex
