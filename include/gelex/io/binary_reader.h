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

#ifndef GELEX_IO_BINARY_READER_H_
#define GELEX_IO_BINARY_READER_H_

#include <cstddef>
#include <cstdint>
#include <filesystem>
#include <format>
#include <string>
#include <string_view>
#include <system_error>
#include <type_traits>
#include <unordered_map>

#include <mio.h>

#include <Eigen/Core>

#include "gelex/exception.h"
#include "gelex/io/binary_format.h"

namespace gelex::detail
{

class BinaryReader
{
   public:
    explicit BinaryReader(std::string_view file_path);

    BinaryReader(const BinaryReader&) = delete;
    BinaryReader(BinaryReader&&) noexcept = default;
    auto operator=(const BinaryReader&) -> BinaryReader& = delete;
    auto operator=(BinaryReader&&) noexcept -> BinaryReader& = default;
    ~BinaryReader() = default;

    [[nodiscard]] auto has_section(
        EffectType effect,
        DataKind kind,
        uint8_t index = 0) const -> bool;

    template <typename eT>
        requires std::is_arithmetic_v<eT>
    [[nodiscard]] auto map(EffectType effect, DataKind kind, uint8_t index = 0)
        const -> Eigen::Map<
            const Eigen::Matrix<eT, Eigen::Dynamic, Eigen::Dynamic>,
            Eigen::Unaligned>;

    template <typename eT>
        requires std::is_arithmetic_v<eT>
    [[nodiscard]] auto mat(EffectType effect, DataKind kind, uint8_t index = 0)
        const -> Eigen::Matrix<eT, Eigen::Dynamic, Eigen::Dynamic>;

    [[nodiscard]] auto n_sections() const -> uint32_t
    {
        return static_cast<uint32_t>(toc_.size());
    }

   private:
    auto parse_footer_and_toc() -> void;

    [[nodiscard]] auto find_entry(
        EffectType effect,
        DataKind kind,
        uint8_t index = 0) const -> const TocEntry&;

    std::filesystem::path path_;
    mio::mmap_source mmap_;
    std::unordered_map<SectionKey, TocEntry, SectionKeyHash> toc_;
};

// --- Implementation ---

inline BinaryReader::BinaryReader(std::string_view file_path)
    : path_(std::string(file_path))
{
    std::error_code ec;
    mmap_.map(path_.string(), ec);
    if (ec)
    {
        if (!std::filesystem::exists(path_))
        {
            throw FileNotFoundException(
                std::format("{}: not found", path_.string()));
        }
        throw FileOpenException(
            std::format(
                "{}: failed to mmap: {}", path_.string(), ec.message()));
    }

    if (mmap_.size() < binary_format::kFooterSize)
    {
        throw FileFormatException(
            std::format(
                "{}: file too small for container footer", path_.string()));
    }

    parse_footer_and_toc();
}

inline auto BinaryReader::parse_footer_and_toc() -> void
{
    const auto file_size = mmap_.size();
    const auto* data = reinterpret_cast<const std::byte*>(mmap_.data());
    const auto path_str = path_.string();

    // Read footer (last 32 bytes)
    const auto* footer = data + file_size - binary_format::kFooterSize;

    binary_format::validate_magic(
        footer, binary_format::kBinaryFormatMagic, path_str, "container");

    const auto version = binary_format::decode<uint32_t>(footer + 8);
    if (version != binary_format::kBinaryFormatVersion)
    {
        throw FileFormatException(
            std::format(
                "{}: unsupported container version {}, expected {}",
                path_str,
                version,
                binary_format::kBinaryFormatVersion));
    }

    const auto toc_offset = binary_format::decode<uint64_t>(footer + 12);
    const auto n_entries = binary_format::decode<uint32_t>(footer + 20);

    // Validate TOC region fits in file
    const auto toc_region_size
        = static_cast<uint64_t>(n_entries) * binary_format::kTocEntrySize;
    if (toc_offset + toc_region_size + binary_format::kFooterSize != file_size)
    {
        throw FileFormatException(
            std::format("{}: TOC region does not match file size", path_str));
    }

    // Read TOC entries
    const auto* toc_data = data + toc_offset;

    for (uint32_t i = 0; i < n_entries; ++i)
    {
        const auto* entry_buf
            = toc_data + static_cast<size_t>(i) * binary_format::kTocEntrySize;

        TocEntry entry{
            .key
            = SectionKey{.effect = static_cast<EffectType>(entry_buf[0]), .kind = static_cast<DataKind>(entry_buf[1]), .index = static_cast<uint8_t>(entry_buf[3])},
            .dtype = static_cast<uint8_t>(entry_buf[2]),
            .offset = binary_format::decode<uint64_t>(entry_buf + 8),
            .size = binary_format::decode<uint64_t>(entry_buf + 16),
            .rows = binary_format::decode<uint64_t>(entry_buf + 24),
            .cols = binary_format::decode<uint64_t>(entry_buf + 32)};

        // Validate section data fits in file
        if (entry.offset + entry.size > toc_offset)
        {
            throw FileFormatException(
                std::format(
                    "{}: section {} data exceeds TOC boundary", path_str, i));
        }

        toc_.emplace(entry.key, entry);
    }
}

inline auto BinaryReader::has_section(
    EffectType effect,
    DataKind kind,
    uint8_t index) const -> bool
{
    return toc_.contains(
        SectionKey{.effect = effect, .kind = kind, .index = index});
}

inline auto BinaryReader::find_entry(
    EffectType effect,
    DataKind kind,
    uint8_t index) const -> const TocEntry&
{
    const SectionKey key{.effect = effect, .kind = kind, .index = index};

    auto it = toc_.find(key);
    if (it == toc_.end())
    {
        throw FileFormatException(
            std::format(
                "{}: section not found (effect={}, kind={})",
                path_.string(),
                static_cast<int>(effect),
                static_cast<int>(kind)));
    }
    return it->second;
}

template <typename eT>
    requires std::is_arithmetic_v<eT>
auto BinaryReader::map(EffectType effect, DataKind kind, uint8_t index) const
    -> Eigen::Map<
        const Eigen::Matrix<eT, Eigen::Dynamic, Eigen::Dynamic>,
        Eigen::Unaligned>
{
    const auto& entry = find_entry(effect, kind, index);
    const auto path_str = path_.string();

    if (entry.dtype != binary_format::dtype_code<eT>())
    {
        throw ArgumentValidationException(
            std::format(
                "{}: dtype mismatch, section={}, requested={}",
                path_str,
                entry.dtype,
                binary_format::dtype_code<eT>()));
    }

    const auto expected_bytes = binary_format::checked_mul(
        binary_format::checked_mul(
            static_cast<size_t>(entry.rows),
            static_cast<size_t>(entry.cols),
            path_str,
            "number of elements"),
        sizeof(eT),
        path_str,
        "payload bytes");

    if (expected_bytes != entry.size)
    {
        throw FileFormatException(
            std::format(
                "{}: section payload size mismatch, expected {} bytes, "
                "got {}",
                path_str,
                expected_bytes,
                entry.size));
    }

    const auto* data_ptr
        = reinterpret_cast<const eT*>(mmap_.data() + entry.offset);

    return Eigen::Map<
        const Eigen::Matrix<eT, Eigen::Dynamic, Eigen::Dynamic>,
        Eigen::Unaligned>(
        data_ptr,
        static_cast<Eigen::Index>(entry.rows),
        static_cast<Eigen::Index>(entry.cols));
}

template <typename eT>
    requires std::is_arithmetic_v<eT>
auto BinaryReader::mat(EffectType effect, DataKind kind, uint8_t index) const
    -> Eigen::Matrix<eT, Eigen::Dynamic, Eigen::Dynamic>
{
    return Eigen::Matrix<eT, Eigen::Dynamic, Eigen::Dynamic>(
        map<eT>(effect, kind, index));
}

}  // namespace gelex::detail

#endif  // GELEX_IO_BINARY_READER_H_
