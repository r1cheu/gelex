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

#include "mmap_reader.h"

#include <format>
#include <system_error>

#include "gelex/exception.h"

namespace gelex::detail
{

namespace
{

constexpr size_t kBedHeaderSize = 3;
constexpr uint8_t kBedMagic0 = 0x6C;
constexpr uint8_t kBedMagic1 = 0x1B;
constexpr uint8_t kBedMagic2 = 0x01;

auto has_valid_bed_magic(const mio::mmap_source& mmap) -> bool
{
    return mmap[0] == kBedMagic0 && mmap[1] == kBedMagic1
           && mmap[2] == kBedMagic2;
}

}  // namespace

BedMmapReader::BedMmapReader(
    const std::filesystem::path& bed_path,
    Eigen::Index num_raw_snps,
    Eigen::Index bytes_per_variant)
    : bytes_per_variant_(bytes_per_variant)
{
    if (bytes_per_variant_ <= 0)
    {
        throw FileFormatException(
            std::format("{}: invalid bytes per variant", bed_path.string()));
    }

    std::error_code ec;
    mmap_.map(bed_path.string(), ec);
    if (ec)
    {
        throw FileOpenException(
            std::format("{}: failed to mmap bed file", bed_path.string()));
    }

    if (mmap_.size() <= kBedHeaderSize)
    {
        throw FileFormatException(
            std::format("{}: bed file too short", bed_path.string()));
    }

    if (!has_valid_bed_magic(mmap_))
    {
        throw FileFormatException(
            std::format("{}: invalid BED magic number", bed_path.string()));
    }

    const size_t expected_size = kBedHeaderSize
                                 + (static_cast<size_t>(num_raw_snps)
                                    * static_cast<size_t>(bytes_per_variant_));
    if (mmap_.size() < expected_size)
    {
        throw FileFormatException(
            std::format(
                "{}: bed file truncated. Expected {} bytes, got {}",
                bed_path.string(),
                expected_size,
                mmap_.size()));
    }

    payload_ptr_
        = reinterpret_cast<const uint8_t*>(mmap_.data()) + kBedHeaderSize;

    const size_t payload_size_bytes = mmap_.size() - kBedHeaderSize;
    payload_num_snps_
        = payload_size_bytes / static_cast<size_t>(bytes_per_variant_);
}

auto BedMmapReader::chunk_ptr(Eigen::Index start_snp, Eigen::Index num_snps)
    const -> const uint8_t*
{
    if (start_snp < 0 || num_snps < 0)
    {
        return nullptr;
    }

    const auto start = static_cast<size_t>(start_snp);
    const auto count = static_cast<size_t>(num_snps);

    if (start > payload_num_snps_)
    {
        return nullptr;
    }

    if (count > (payload_num_snps_ - start))
    {
        return nullptr;
    }

    const auto offset = start * static_cast<size_t>(bytes_per_variant_);
    return payload_ptr_ + offset;
}

}  // namespace gelex::detail
