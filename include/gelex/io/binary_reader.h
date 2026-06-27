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

#include <fmt/format.h>
#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <filesystem>
#include <string>
#include <string_view>
#include <system_error>
#include <type_traits>
#include <unordered_map>
#include <vector>

#include <Eigen/Core>

#include "gelex/exception.h"
#include "gelex/infra/string_hash.h"
#include "gelex/io/detail/binary_format.h"
#include "gelex/io/mapped_file.h"

namespace gelex::io
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

    [[nodiscard]] auto contains(std::string_view path) const -> bool;

    template <typename eT>
        requires std::is_arithmetic_v<eT>
    [[nodiscard]] auto to_map(std::string_view path) const -> Eigen::Map<
        const Eigen::Matrix<eT, Eigen::Dynamic, Eigen::Dynamic>,
        Eigen::Unaligned>;

    template <typename eT>
        requires std::is_arithmetic_v<eT>
    [[nodiscard]] auto to_mat(std::string_view path) const
        -> Eigen::Matrix<eT, Eigen::Dynamic, Eigen::Dynamic>;

    template <typename eT>
        requires std::is_arithmetic_v<eT>
    [[nodiscard]] auto is_type(std::string_view path) const -> bool
    {
        return find_entry(path).dtype == detail::TYPE_BYTE<eT>;
    }

    [[nodiscard]] auto to_strings(std::string_view path) const
        -> std::vector<std::string_view>;

    [[nodiscard]] auto n_sections() const -> uint64_t
    {
        return static_cast<uint64_t>(toc_.size());
    }

    [[nodiscard]] auto section_paths() const -> std::vector<std::string_view>
    {
        std::vector<std::string_view> paths;
        paths.reserve(toc_.size());
        for (const auto& [key, _] : toc_)
        {
            paths.emplace_back(key);
        }
        return paths;
    }

   private:
    auto parse_footer_and_toc() -> void;

    [[nodiscard]] auto find_entry(std::string_view path) const
        -> const detail::TocEntry&;

    std::filesystem::path path_;
    MappedFile mmap_;
    std::unordered_map<
        std::string,
        detail::TocEntry,
        infra::TransparentHash<std::string>,
        infra::TransparentEqual<std::string>>
        toc_;
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
            throw GelexException(fmt::format("{}: not found", path_.string()));
        }
        throw GelexException(
            fmt::format(
                "{}: failed to mmap: {}", path_.string(), ec.message()));
    }

    if (mmap_.size() < detail::FOOTER_SIZE)
    {
        throw GelexException(
            fmt::format(
                "{}: file too small for container footer", path_.string()));
    }

    parse_footer_and_toc();
}

inline auto BinaryReader::parse_footer_and_toc() -> void
{
    const auto file_size = mmap_.size();
    const auto* data = mmap_.data();
    const auto path_str = path_.string();

    const auto* footer = data + file_size - detail::FOOTER_SIZE;

    if (!std::equal(
            detail::BINARY_FORMAT_MAGIC.begin(),
            detail::BINARY_FORMAT_MAGIC.end(),
            footer))
    {
        throw GelexException(
            fmt::format("{}: invalid container magic", path_str));
    }

    const auto toc_offset = detail::decode<uint64_t>(footer + 8);
    const auto n_entries = detail::decode<uint64_t>(footer + 16);

    const auto toc_region_size
        = static_cast<uint64_t>(n_entries) * detail::TOC_ENTRY_SIZE;
    if (toc_offset + toc_region_size + detail::FOOTER_SIZE != file_size)
    {
        throw GelexException(
            fmt::format("{}: TOC region does not match file size", path_str));
    }

    const auto* toc_data = data + toc_offset;

    for (uint64_t i = 0; i < n_entries; ++i)
    {
        const auto* entry_buf
            = toc_data + static_cast<size_t>(i) * detail::TOC_ENTRY_SIZE;
        auto entry = detail::TocEntry::from_bytes(entry_buf);

        if (entry.offset + entry.size > toc_offset)
        {
            throw GelexException(
                fmt::format(
                    "{}: section {} data exceeds TOC boundary", path_str, i));
        }

        auto key = std::string(detail::path_as_view(entry.path));
        toc_.emplace(std::move(key), entry);
    }
}

inline auto BinaryReader::contains(std::string_view path) const -> bool
{
    return toc_.contains(path);
}

inline auto BinaryReader::find_entry(std::string_view path) const
    -> const detail::TocEntry&
{
    auto it = toc_.find(path);
    if (it == toc_.end())
    {
        throw GelexException(
            fmt::format("{}: section not found: \"{}\"", path_.string(), path));
    }
    return it->second;
}

inline auto BinaryReader::to_strings(std::string_view path) const
    -> std::vector<std::string_view>
{
    const auto& entry = find_entry(path);

    if (entry.dtype != detail::TYPE_STRING)
    {
        throw GelexException(
            fmt::format(
                "{}: section \"{}\" is not a string section (dtype={})",
                path_.string(),
                path,
                entry.dtype));
    }

    const auto* base
        = reinterpret_cast<const char*>(mmap_.data()) + entry.offset;
    const auto* end = base + entry.size;

    std::vector<std::string_view> result;
    result.reserve(entry.rows);

    const auto* cursor = base;
    while (cursor < end)
    {
        const auto* null_pos = reinterpret_cast<const char*>(
            std::memchr(cursor, '\0', static_cast<size_t>(end - cursor)));
        if (null_pos == nullptr)
        {
            throw GelexException(
                fmt::format(
                    "{}: string section \"{}\" missing null terminator",
                    path_.string(),
                    path));
        }
        result.emplace_back(cursor, static_cast<size_t>(null_pos - cursor));
        cursor = null_pos + 1;
    }

    if (result.size() != entry.rows)
    {
        throw GelexException(
            fmt::format(
                "{}: string section \"{}\" count mismatch: expected {}, "
                "got {}",
                path_.string(),
                path,
                entry.rows,
                result.size()));
    }

    return result;
}

template <typename eT>
    requires std::is_arithmetic_v<eT>
auto BinaryReader::to_map(std::string_view path) const -> Eigen::Map<
    const Eigen::Matrix<eT, Eigen::Dynamic, Eigen::Dynamic>,
    Eigen::Unaligned>
{
    const auto& entry = find_entry(path);
    const auto path_str = path_.string();

    if (entry.dtype != detail::TYPE_BYTE<eT>)
    {
        throw GelexException(
            fmt::format(
                "{}: dtype mismatch, section={}, requested={}",
                path_str,
                entry.dtype,
                detail::TYPE_BYTE<eT>));
    }

    const auto n_elements
        = static_cast<size_t>(entry.rows) * static_cast<size_t>(entry.cols);
    const auto expected_bytes = n_elements * sizeof(eT);

    if (expected_bytes != entry.size)
    {
        throw GelexException(
            fmt::format(
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
auto BinaryReader::to_mat(std::string_view path) const
    -> Eigen::Matrix<eT, Eigen::Dynamic, Eigen::Dynamic>
{
    const auto& entry = find_entry(path);
    if (entry.dtype == detail::TYPE_BYTE<eT>)
    {
        return Eigen::Matrix<eT, Eigen::Dynamic, Eigen::Dynamic>(
            to_map<eT>(path));
    }

    auto cast_from = [&]<typename StoredT>()
        -> Eigen::Matrix<eT, Eigen::Dynamic, Eigen::Dynamic>
    { return to_map<StoredT>(path).template cast<eT>(); };

    switch (entry.dtype)
    {
        case detail::TYPE_BYTE<double>:
            return cast_from.template operator()<double>();
        case detail::TYPE_BYTE<int8_t>:
            return cast_from.template operator()<int8_t>();
        case detail::TYPE_BYTE<uint8_t>:
            return cast_from.template operator()<uint8_t>();
        case detail::TYPE_BYTE<int64_t>:
            return cast_from.template operator()<int64_t>();
        default:
            throw GelexException(
                fmt::format(
                    "{}: unsupported dtype cast, section={}, dtype={}",
                    path_.string(),
                    path,
                    entry.dtype));
    }
}

}  // namespace gelex::io

#endif  // GELEX_IO_BINARY_READER_H_
