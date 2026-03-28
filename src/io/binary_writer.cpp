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

#include "gelex/io/binary_writer.h"

#include <algorithm>
#include <array>
#include <cstddef>
#include <cstdint>
#include <format>
#include <fstream>
#include <string>
#include <string_view>

#include "gelex/exception.h"
#include "gelex/io/binary_format.h"
#include "gelex/io/parser.h"

namespace gelex::detail
{

BinaryWriter::BinaryWriter(std::string_view output_path)
    : output_path_(std::string(output_path)),
      file_(
          open_file<std::ofstream>(
              output_path_,
              std::ios::binary | std::ios::trunc))
{
}

auto BinaryWriter::check_duplicate_path(std::string_view path) const -> void
{
    for (const auto& rs : reserved_)
    {
        if (binary_format::path_as_view(rs.entry.path) == path)
        {
            throw GelexException(
                std::format(
                    "{}: duplicate section path \"{}\"",
                    output_path_.string(),
                    path));
        }
    }
}

auto BinaryWriter::reserve(
    std::string_view path,
    uint8_t dtype,
    uint64_t rows,
    uint64_t cols) -> size_t
{
    if (path.size() > binary_format::kMaxPathLength)
    {
        throw GelexException(
            std::format(
                "{}: path too long ({} > {}): \"{}\"",
                output_path_.string(),
                path.size(),
                binary_format::kMaxPathLength,
                path));
    }

    check_duplicate_path(path);

    TocEntry entry;
    std::copy(path.begin(), path.end(), entry.path.begin());
    entry.dtype = dtype;

    const auto element_size
        = static_cast<uint64_t>(static_cast<unsigned>(dtype) >> 2U);
    const auto total_bytes = rows * cols * element_size;
    const auto aligned_offset
        = align_up(next_offset_, binary_format::kPageAlignment);
    entry.offset = aligned_offset;
    entry.size = total_bytes;
    entry.rows = rows;
    entry.cols = cols;

    next_offset_ = aligned_offset + total_bytes;
    reserved_.push_back(
        ReservedSection{.entry = entry, .cursor = aligned_offset});

    return reserved_.size() - 1;
}

auto BinaryWriter::write_raw(
    size_t handle,
    const char* data,
    std::streamsize bytes) -> void
{
    if (handle >= reserved_.size())
    {
        throw GelexException(
            std::format(
                "{}: invalid section handle {}",
                output_path_.string(),
                handle));
    }

    auto& rs = reserved_[handle];
    const auto end_bound = rs.entry.offset + rs.entry.size;
    if (rs.cursor + static_cast<uint64_t>(bytes) > end_bound)
    {
        throw GelexException(
            std::format(
                "{}: write overflow: cursor={}, bytes={}, limit={}",
                output_path_.string(),
                rs.cursor,
                bytes,
                end_bound));
    }

    if (rs.cursor != file_cursor_)
    {
        file_.seekp(static_cast<std::streamoff>(rs.cursor));
    }
    file_.write(data, bytes);
    rs.cursor += static_cast<uint64_t>(bytes);
    file_cursor_ = rs.cursor;
}

auto BinaryWriter::align_up(uint64_t value, uint64_t alignment) noexcept
    -> uint64_t
{
    return (value + alignment - 1) / alignment * alignment;
}

auto BinaryWriter::write_footer(uint64_t toc_offset, uint64_t n_sections)
    -> void
{
    std::array<std::byte, binary_format::kFooterSize> buf{};
    std::copy(
        binary_format::kBinaryFormatMagic.begin(),
        binary_format::kBinaryFormatMagic.end(),
        buf.begin());
    binary_format::encode(toc_offset, &buf[8]);
    binary_format::encode(n_sections, &buf[16]);

    file_.write(
        reinterpret_cast<const char*>(buf.data()),
        static_cast<std::streamsize>(buf.size()));
}

auto BinaryWriter::finalize() -> void
{
    for (const auto& rs : reserved_)
    {
        const auto expected = rs.entry.offset + rs.entry.size;
        if (rs.cursor != expected)
        {
            throw GelexException(
                std::format(
                    "{}: section not fully written: cursor={}, "
                    "expected={}",
                    output_path_.string(),
                    rs.cursor,
                    expected));
        }
    }

    const auto toc_offset
        = align_up(next_offset_, binary_format::kPageAlignment);
    file_.seekp(static_cast<std::streamoff>(toc_offset));
    for (const auto& rs : reserved_)
    {
        auto buf = rs.entry.to_bytes();
        file_.write(
            reinterpret_cast<const char*>(buf.data()),
            static_cast<std::streamsize>(buf.size()));
    }

    write_footer(toc_offset, static_cast<uint64_t>(reserved_.size()));

    file_.flush();
    if (!file_.good())
    {
        throw GelexException(
            std::format(
                "{}: failed to write binary file", output_path_.string()));
    }
}

}  // namespace gelex::detail
