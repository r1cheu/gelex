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
#include <exception>
#include <fmt/format.h>
#include <ios>
#include <string>
#include <string_view>

#include "gelex/exception.h"
#include "gelex/infra/log.h"
#include "gelex/io/detail/binary_format.h"

namespace gelex
{

BinaryWriter::BinaryWriter(std::string_view output_path)
    : file_(std::string(output_path), std::ios::binary | std::ios::trunc)
{
}

BinaryWriter::~BinaryWriter() noexcept
{
    if (finalized_)
    {
        return;
    }
    if (std::uncaught_exceptions() > 0)
    {
        return;
    }

    try
    {
        close();
    }
    catch (const std::exception& e)
    {
        gelex::error(
            fmt::format(
                "{}: failed to finalize, discarding output: {}",
                file_.path().string(),
                e.what()));
    }
}

auto BinaryWriter::close() -> void
{
    if (finalized_)
    {
        return;
    }
    finalize();
    file_.commit();
    finalized_ = true;
}

auto BinaryWriter::check_duplicate_path(std::string_view path) const -> void
{
    for (const auto& rs : reserved_)
    {
        if (detail::path_as_view(rs.entry.path) == path)
        {
            throw GelexException(
                fmt::format(
                    "{}: duplicate section path \"{}\"",
                    file_.path().string(),
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
    const auto element_size
        = static_cast<uint64_t>(static_cast<unsigned>(dtype) >> 2U);
    return reserve_section(path, dtype, rows, cols, rows * cols * element_size);
}

auto BinaryWriter::reserve_strings(
    std::string_view path,
    uint64_t rows,
    uint64_t bytes) -> size_t
{
    return reserve_section(path, detail::TYPE_STRING, rows, 0, bytes);
}

auto BinaryWriter::reserve_section(
    std::string_view path,
    uint8_t dtype,
    uint64_t rows,
    uint64_t cols,
    uint64_t bytes) -> size_t
{
    if (path.size() > detail::MAX_PATH_LENGTH)
    {
        throw GelexException(
            fmt::format(
                "{}: path too long ({} > {}): \"{}\"",
                file_.path().string(),
                path.size(),
                detail::MAX_PATH_LENGTH,
                path));
    }

    check_duplicate_path(path);

    detail::TocEntry entry;
    std::copy(path.begin(), path.end(), entry.path.begin());
    entry.dtype = dtype;
    entry.rows = rows;
    entry.cols = cols;
    entry.size = bytes;

    const auto aligned_offset = align_up(next_offset_, detail::PAGE_ALIGNMENT);
    entry.offset = aligned_offset;
    next_offset_ = aligned_offset + bytes;

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
            fmt::format(
                "{}: invalid section handle {}",
                file_.path().string(),
                handle));
    }

    auto& rs = reserved_[handle];
    const auto end_bound = rs.entry.offset + rs.entry.size;
    if (rs.cursor + static_cast<uint64_t>(bytes) > end_bound)
    {
        throw GelexException(
            fmt::format(
                "{}: write overflow: cursor={}, bytes={}, limit={}",
                file_.path().string(),
                rs.cursor,
                bytes,
                end_bound));
    }

    if (rs.cursor != file_cursor_)
    {
        file_.seek(static_cast<std::streamoff>(rs.cursor));
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
    std::array<std::byte, detail::FOOTER_SIZE> buf{};
    std::copy(
        detail::BINARY_FORMAT_MAGIC.begin(),
        detail::BINARY_FORMAT_MAGIC.end(),
        buf.begin());
    detail::encode(toc_offset, &buf[8]);
    detail::encode(n_sections, &buf[16]);

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
                fmt::format(
                    "{}: section not fully written: cursor={}, "
                    "expected={}",
                    file_.path().string(),
                    rs.cursor,
                    expected));
        }
    }

    const auto toc_offset = align_up(next_offset_, detail::PAGE_ALIGNMENT);
    file_.seek(static_cast<std::streamoff>(toc_offset));
    for (const auto& rs : reserved_)
    {
        auto buf = rs.entry.to_bytes();
        file_.write(
            reinterpret_cast<const char*>(buf.data()),
            static_cast<std::streamsize>(buf.size()));
    }

    write_footer(toc_offset, static_cast<uint64_t>(reserved_.size()));
}

}  // namespace gelex
