// Copyright 2026 RuLei Chen
// SPDX-License-Identifier: Apache-2.0

#include "gelex/io/csc_reader.h"

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <filesystem>
#include <fmt/format.h>
#include <ranges>
#include <span>
#include <string>
#include <string_view>
#include <system_error>
#include <utility>
#include <vector>

#include "gelex/exception.h"
#include "gelex/io/binary_format.h"
#include "gelex/io/detail/binary_wire.h"

namespace gelex
{

namespace
{

struct Footer
{
    std::uint64_t index_offset;
    std::uint64_t matrix_count;
};

auto parse_footer(std::span<const std::byte> file, std::string_view path)
    -> Footer
{
    if (file.size() < detail::footer_size)
    {
        throw GelexException(
            fmt::format("{}: file too small for CSC footer", path));
    }
    const auto footer_offset = file.size() - detail::footer_size;
    detail::ByteCursor footer{file.subspan(footer_offset), path};
    if (!std::ranges::equal(
            detail::csc_format_magic,
            footer.read_bytes(detail::csc_format_magic.size(), "magic")))
    {
        throw GelexException(fmt::format("{}: invalid CSC magic", path));
    }
    const auto index_offset = footer.read<std::uint64_t>("footer");
    const auto matrix_count = footer.read<std::uint64_t>("footer");
    if (index_offset % detail::payload_alignment != 0
        || index_offset > footer_offset)
    {
        throw GelexException(fmt::format("{}: invalid CSC index offset", path));
    }
    return Footer{.index_offset = index_offset, .matrix_count = matrix_count};
}

// Each array must be 64-aligned and lie inside [0, index_offset).
auto parse_index_entry(
    detail::ByteCursor& cursor,
    const Footer& footer,
    std::uint64_t number,
    std::string_view path) -> detail::CscEntry
{
    detail::CscEntry entry;
    entry.header = detail::read_matrix_header(cursor);
    entry.nnz = cursor.read<std::uint64_t>("matrix nnz");
    const auto sizes = detail::csc_array_sizes(
        entry.header.type, entry.header.shape[1], entry.nnz);
    for (auto [offset, size] : std::views::zip(entry.offsets, sizes))
    {
        offset = cursor.read<std::uint64_t>("array offset");
        if (offset % detail::payload_alignment != 0
            || offset > footer.index_offset
            || size > footer.index_offset - offset)
        {
            throw GelexException(
                fmt::format(
                    "{}: matrix {} is outside the data region", path, number));
        }
    }
    return entry;
}

}  // namespace

CscReader::CscReader(std::string_view file_path) : path_(std::string(file_path))
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
    parse_footer_and_index();
}

auto CscReader::parse_footer_and_index() -> void
{
    const std::span<const std::byte> file{mmap_.data(), mmap_.size()};
    const auto path = path_.string();
    const auto footer = parse_footer(file, path);
    const auto footer_offset = file.size() - detail::footer_size;
    const auto index_offset = static_cast<std::size_t>(footer.index_offset);
    detail::ByteCursor index{
        file.subspan(index_offset, footer_offset - index_offset), path};
    if (footer.matrix_count > index.size() / detail::csc_index_entry_fixed_size)
    {
        throw GelexException(fmt::format("{}: invalid CSC matrix count", path));
    }

    std::vector<std::pair<std::uint64_t, std::uint64_t>> ranges;
    for (std::uint64_t number = 0; number < footer.matrix_count; ++number)
    {
        auto entry = parse_index_entry(index, footer, number, path);
        const auto sizes = detail::csc_array_sizes(
            entry.header.type, entry.header.shape[1], entry.nnz);
        for (auto [offset, size] : std::views::zip(entry.offsets, sizes))
        {
            ranges.emplace_back(offset, offset + size);
        }
        entries_.push_back(std::move(entry));
    }
    if (!index.empty())
    {
        throw GelexException(fmt::format("{}: CSC index size mismatch", path));
    }

    index_.reserve(entries_.size());
    for (const auto& entry : entries_)
    {
        const auto [_, inserted]
            = index_.emplace(entry.header.identifier, index_.size());
        if (!inserted)
        {
            throw GelexException(
                fmt::format(
                    "{}: duplicate matrix identifier \"{}\"",
                    path,
                    entry.header.identifier));
        }
    }

    std::ranges::sort(ranges);
    for (std::size_t number = 1; number < ranges.size(); ++number)
    {
        if (ranges[number - 1].second > ranges[number].first)
        {
            throw GelexException(fmt::format("{}: CSC arrays overlap", path));
        }
    }
}

auto CscReader::contains(std::string_view identifier) const -> bool
{
    return index_.contains(identifier);
}

auto CscReader::size() const noexcept -> std::size_t
{
    return entries_.size();
}

auto CscReader::info(std::string_view identifier) const -> const MatrixHeader&
{
    return find_entry(identifier).header;
}

auto CscReader::nnz(std::string_view identifier) const -> std::uint64_t
{
    return find_entry(identifier).nnz;
}

auto CscReader::payloads() const -> std::vector<MatrixHeader>
{
    std::vector<MatrixHeader> result;
    result.reserve(entries_.size());
    for (const auto& entry : entries_)
    {
        result.push_back(entry.header);
    }
    std::ranges::sort(result, {}, &MatrixHeader::identifier);
    return result;
}

auto CscReader::find_entry(std::string_view identifier) const
    -> const detail::CscEntry&
{
    const auto iterator = index_.find(identifier);
    if (iterator == index_.end())
    {
        throw GelexException(
            fmt::format(
                "{}: matrix not found: \"{}\"", path_.string(), identifier));
    }
    return entries_[iterator->second];
}

}  // namespace gelex
