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

#include <algorithm>
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
        gelex::EffectType effect,
        DataKind kind,
        uint8_t index = 0) const -> bool;

    template <typename eT>
        requires std::is_arithmetic_v<eT>
    [[nodiscard]] auto
    map(gelex::EffectType effect, DataKind kind, uint8_t index = 0) const
        -> Eigen::Map<
            const Eigen::Matrix<eT, Eigen::Dynamic, Eigen::Dynamic>,
            Eigen::Unaligned>;

    template <typename eT>
        requires std::is_arithmetic_v<eT>
    [[nodiscard]] auto
    mat(gelex::EffectType effect, DataKind kind, uint8_t index = 0) const
        -> Eigen::Matrix<eT, Eigen::Dynamic, Eigen::Dynamic>;

    [[nodiscard]] auto n_sections() const -> uint64_t
    {
        return static_cast<uint64_t>(toc_.size());
    }

   private:
    auto parse_footer_and_toc() -> void;

    [[nodiscard]] auto find_entry(
        gelex::EffectType effect,
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

    if (!std::equal(
            binary_format::kBinaryFormatMagic.begin(),
            binary_format::kBinaryFormatMagic.end(),
            footer))
    {
        throw FileFormatException(
            std::format("{}: invalid container magic", path_str));
    }

    const auto toc_offset = binary_format::decode<uint64_t>(footer + 8);
    const auto n_entries = binary_format::decode<uint64_t>(footer + 16);

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

    for (uint64_t i = 0; i < n_entries; ++i)
    {
        const auto* entry_buf
            = toc_data + static_cast<size_t>(i) * binary_format::kTocEntrySize;
        auto entry = TocEntry::from_bytes(entry_buf);

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
    gelex::EffectType effect,
    DataKind kind,
    uint8_t index) const -> bool
{
    return toc_.contains(
        SectionKey{.effect = effect, .kind = kind, .index = index});
}

inline auto BinaryReader::find_entry(
    gelex::EffectType effect,
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
                static_cast<int>(effect.to_byte()),
                static_cast<int>(kind)));
    }
    return it->second;
}

template <typename eT>
    requires std::is_arithmetic_v<eT>
auto BinaryReader::map(gelex::EffectType effect, DataKind kind, uint8_t index)
    const -> Eigen::Map<
        const Eigen::Matrix<eT, Eigen::Dynamic, Eigen::Dynamic>,
        Eigen::Unaligned>
{
    const auto& entry = find_entry(effect, kind, index);
    const auto path_str = path_.string();

    if (entry.dtype != binary_format::kTypeByte<eT>)
    {
        throw ArgumentValidationException(
            std::format(
                "{}: dtype mismatch, section={}, requested={}",
                path_str,
                entry.dtype,
                binary_format::kTypeByte<eT>));
    }

    const auto n_elements
        = static_cast<size_t>(entry.rows) * static_cast<size_t>(entry.cols);
    const auto expected_bytes = n_elements * sizeof(eT);

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
auto BinaryReader::mat(gelex::EffectType effect, DataKind kind, uint8_t index)
    const -> Eigen::Matrix<eT, Eigen::Dynamic, Eigen::Dynamic>
{
    return Eigen::Matrix<eT, Eigen::Dynamic, Eigen::Dynamic>(
        map<eT>(effect, kind, index));
}

}  // namespace gelex::detail

#endif  // GELEX_IO_BINARY_READER_H_
