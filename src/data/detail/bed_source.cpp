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

#include "gelex/data/detail/bed_source.h"

#include <cstddef>
#include <cstdint>
#include <filesystem>
#include <fmt/format.h>
#include <string>
#include <system_error>

#include "gelex/io/detail/parser.h"

namespace
{

inline constexpr std::size_t bed_header_size = 3;
inline constexpr std::uint8_t bed_magic_0 = 0x6C;
inline constexpr std::uint8_t bed_magic_1 = 0x1B;
inline constexpr std::uint8_t bed_magic_2 = 0x01;

}  // namespace

namespace gelex::detail
{

BedSource::BedSource(
    const std::filesystem::path& bed_path,
    size_type num_variants,
    size_type num_samples)
    : num_variants_(num_variants),
      num_samples_(num_samples),
      bytes_per_variant_((num_samples + 3) / 4)
{
    if (num_variants == 0)
    {
        throw GelexException(
            fmt::format("{}: BED variant count is zero", bed_path.string()));
    }

    if (num_samples == 0)
    {
        throw GelexException(
            fmt::format("{}: BED sample count is zero", bed_path.string()));
    }

    std::error_code ec;
    mmap_.map(bed_path.string(), ec);
    if (ec)
    {
        throw GelexException(
            fmt::format("{}: failed to mmap bed file", bed_path.string()));
    }

    if (mmap_.size() < bed_header_size)
    {
        throw GelexException(
            fmt::format("{}: bed file too short", bed_path.string()));
    }

    const auto* raw = reinterpret_cast<const std::uint8_t*>(mmap_.data());
    if (raw[0] != bed_magic_0 || raw[1] != bed_magic_1 || raw[2] != bed_magic_2)
    {
        throw GelexException(
            fmt::format("{}: invalid BED magic number", bed_path.string()));
    }

    const auto expected_size
        = bed_header_size + (num_variants_ * bytes_per_variant_);
    if (mmap_.size() != expected_size)
    {
        throw GelexException(
            fmt::format(
                "{}: bed file size mismatch. Expected {} bytes, got {}",
                bed_path.string(),
                expected_size,
                mmap_.size()));
    }

    payload_ = raw + bed_header_size;
}

auto open_bed_source(const std::string& bfile_prefix) -> BedSource
{
    const std::filesystem::path fam_path = bfile_prefix + ".fam";
    const std::filesystem::path bim_path = bfile_prefix + ".bim";
    const std::filesystem::path bed_path = bfile_prefix + ".bed";

    const auto num_samples = detail::count_total_lines(fam_path);
    const auto num_variants = detail::count_total_lines(bim_path);

    return BedSource{bed_path, num_variants, num_samples};
}

}  // namespace gelex::detail
