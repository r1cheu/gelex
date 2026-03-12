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

#include <array>
#include <cstddef>
#include <cstdint>
#include <format>
#include <fstream>
#include <string>

#include "gelex/exception.h"
#include "gelex/io/binary_format.h"

namespace gelex::detail
{

namespace
{

auto align_up(uint64_t value, uint64_t alignment) -> uint64_t
{
    return (value + alignment - 1) / alignment * alignment;
}

auto write_toc_entry(std::ofstream& out, const TocEntry& entry) -> void
{
    std::array<std::byte, binary_format::kTocEntrySize> buf{};
    buf[0] = static_cast<std::byte>(entry.key.effect);
    buf[1] = static_cast<std::byte>(entry.key.kind);
    buf[2] = static_cast<std::byte>(entry.dtype);
    buf[3] = static_cast<std::byte>(entry.key.index);
    // buf[4..7] = padding (already zero)
    binary_format::encode(entry.offset, &buf[8]);
    binary_format::encode(entry.size, &buf[16]);
    binary_format::encode(entry.rows, &buf[24]);
    binary_format::encode(entry.cols, &buf[32]);

    out.write(
        reinterpret_cast<const char*>(buf.data()),
        static_cast<std::streamsize>(buf.size()));
}

auto write_footer(std::ofstream& out, uint64_t toc_offset, uint32_t n_sections)
    -> void
{
    std::array<std::byte, binary_format::kFooterSize> buf{};
    std::copy(
        binary_format::kBinaryFormatMagic.begin(),
        binary_format::kBinaryFormatMagic.end(),
        buf.begin());
    binary_format::encode(binary_format::kBinaryFormatVersion, &buf[8]);
    binary_format::encode(toc_offset, &buf[12]);
    binary_format::encode(n_sections, &buf[20]);
    // buf[24..31] = reserved padding (already zero)

    out.write(
        reinterpret_cast<const char*>(buf.data()),
        static_cast<std::streamsize>(buf.size()));
}

}  // namespace

BinaryWriter::BinaryWriter(std::string_view output_path)
    : output_path_(std::string(output_path))
{
}

auto BinaryWriter::check_duplicate_key(const SectionKey& key) -> void
{
    for (const auto& rs : reserved_)
    {
        if (rs.entry.key == key)
        {
            throw ArgumentValidationException(
                std::format(
                    "{}: duplicate section key (effect={}, kind={}, index={})",
                    output_path_.string(),
                    static_cast<int>(key.effect),
                    static_cast<int>(key.kind),
                    key.index));
        }
    }
}

auto BinaryWriter::open() -> std::ofstream&
{
    if (!file_.is_open())
    {
        file_.open(output_path_, std::ios::binary | std::ios::trunc);
        if (!file_.is_open())
        {
            throw FileOpenException(
                std::format(
                    "{}: failed to open container output",
                    output_path_.string()));
        }
    }
    return file_;
}

auto BinaryWriter::reserve(
    SectionKey key,
    uint8_t dtype,
    uint64_t rows,
    uint64_t cols) -> size_t
{
    check_duplicate_key(key);
    open();

    const auto element_size
        = static_cast<uint64_t>(static_cast<unsigned>(dtype) >> 2U);
    const auto total_bytes = rows * cols * element_size;
    const auto aligned_offset
        = align_up(next_offset_, binary_format::kPageAlignment);
    next_offset_ = aligned_offset + total_bytes;

    reserved_.push_back(
        ReservedSection{
            .entry
            = TocEntry{.key = key, .dtype = dtype, .offset = aligned_offset, .size = total_bytes, .rows = rows, .cols = cols},
            .cursor = aligned_offset});

    return reserved_.size() - 1;
}

auto BinaryWriter::write(size_t handle, const char* data, std::streamsize bytes)
    -> void
{
    if (handle >= reserved_.size())
    {
        throw ArgumentValidationException(
            std::format(
                "{}: invalid section handle {}",
                output_path_.string(),
                handle));
    }

    auto& rs = reserved_[handle];
    const auto end_bound = rs.entry.offset + rs.entry.size;
    if (rs.cursor + static_cast<uint64_t>(bytes) > end_bound)
    {
        throw ArgumentValidationException(
            std::format(
                "{}: write overflow: cursor={}, bytes={}, limit={}",
                output_path_.string(),
                rs.cursor,
                bytes,
                end_bound));
    }

    file_.seekp(static_cast<std::streamoff>(rs.cursor));
    file_.write(data, bytes);
    rs.cursor += static_cast<uint64_t>(bytes);
}

auto BinaryWriter::finalize() -> void
{
    for (const auto& rs : reserved_)
    {
        const auto expected = rs.entry.offset + rs.entry.size;
        if (rs.cursor != expected)
        {
            throw FileWriteException(
                std::format(
                    "{}: section not fully written: cursor={}, "
                    "expected={}",
                    output_path_.string(),
                    rs.cursor,
                    expected));
        }
    }

    auto& out = open();

    const auto toc_pos = align_up(next_offset_, binary_format::kPageAlignment);
    out.seekp(static_cast<std::streamoff>(toc_pos));

    const auto toc_offset = static_cast<uint64_t>(out.tellp());
    for (const auto& rs : reserved_)
    {
        write_toc_entry(out, rs.entry);
    }

    write_footer(out, toc_offset, static_cast<uint32_t>(reserved_.size()));

    out.flush();
    if (!out.good())
    {
        throw FileWriteException(
            std::format(
                "{}: failed to write container", output_path_.string()));
    }
}

}  // namespace gelex::detail
