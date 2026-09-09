// Copyright 2026 RuLei Chen
// SPDX-License-Identifier: Apache-2.0

#include "gelex/io/dense_reader.h"

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <filesystem>
#include <fmt/format.h>
#include <span>
#include <string>
#include <string_view>
#include <system_error>
#include <utility>
#include <vector>

#include "gelex/exception.h"
#include "gelex/io/detail/binary_wire.h"

namespace gelex
{

namespace
{

struct Footer
{
    std::uint64_t directory_offset;
    std::uint64_t payload_count;
};

auto parse_footer(std::span<const std::byte> file, std::string_view path)
    -> Footer
{
    if (file.size() < detail::footer_size)
    {
        throw GelexException(
            fmt::format("{}: file too small for container footer", path));
    }

    const auto footer_offset = file.size() - detail::footer_size;
    detail::ByteCursor footer{file.subspan(footer_offset), path};
    if (!std::ranges::equal(
            detail::binary_format_magic,
            footer.read_bytes(detail::binary_format_magic.size(), "magic")))
    {
        throw GelexException(fmt::format("{}: invalid container magic", path));
    }
    const auto directory_offset = footer.read<std::uint64_t>("footer");
    const auto payload_count = footer.read<std::uint64_t>("footer");
    if (directory_offset % detail::payload_alignment != 0
        || directory_offset > footer_offset)
    {
        throw GelexException(fmt::format("{}: invalid directory offset", path));
    }

    return Footer{
        .directory_offset = directory_offset, .payload_count = payload_count};
}

auto parse_payload_entry(
    detail::ByteCursor& cursor,
    const Footer& footer,
    std::uint64_t index,
    std::string_view path) -> detail::MatrixEntry
{
    auto header = detail::read_matrix_header(cursor);
    const auto offset = cursor.read<std::uint64_t>("payload offset");
    const auto size = cursor.read<std::uint64_t>("payload size");

    if (offset % detail::payload_alignment != 0
        || offset > footer.directory_offset
        || size > footer.directory_offset - offset)
    {
        throw GelexException(
            fmt::format(
                "{}: payload {} is outside the data region", path, index));
    }

    const auto element_count
        = detail::checked_product(header.shape[0], header.shape[1]);
    if (size
        != detail::checked_product(
            element_count, detail::binary_type_size(header.type)))
    {
        throw GelexException(
            fmt::format(
                "{}: payload {} size does not match shape", path, index));
    }

    return detail::MatrixEntry{
        .header = std::move(header), .offset = offset, .size = size};
}

}  // namespace

DenseReader::DenseReader(std::string_view file_path)
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

    parse_footer_and_index();
}

auto DenseReader::parse_footer_and_index() -> void
{
    const std::span<const std::byte> file{mmap_.data(), mmap_.size()};
    const auto path_string = path_.string();
    const auto footer = parse_footer(file, path_string);
    const auto footer_offset = file.size() - detail::footer_size;
    const auto directory_offset
        = static_cast<std::size_t>(footer.directory_offset);
    detail::ByteCursor directory{
        file.subspan(directory_offset, footer_offset - directory_offset),
        path_string};
    if (footer.payload_count
        > directory.size() / detail::payload_entry_fixed_size)
    {
        throw GelexException(
            fmt::format("{}: invalid payload count", path_string));
    }

    std::vector<std::pair<std::uint64_t, std::uint64_t>> payload_ranges;

    for (std::uint64_t index = 0; index < footer.payload_count; ++index)
    {
        auto entry = parse_payload_entry(directory, footer, index, path_string);
        if (entry.size != 0)
        {
            payload_ranges.emplace_back(
                entry.offset, entry.offset + entry.size);
        }
        payloads_.push_back(std::move(entry));
    }

    if (!directory.empty())
    {
        throw GelexException(
            fmt::format("{}: directory size mismatch", path_string));
    }

    index_.reserve(payloads_.size());
    for (const auto& payload : payloads_)
    {
        const auto [_, inserted]
            = index_.emplace(payload.header.identifier, index_.size());
        if (!inserted)
        {
            throw GelexException(
                fmt::format(
                    "{}: duplicate payload identifier \"{}\"",
                    path_string,
                    payload.header.identifier));
        }
    }

    validate_payload_ranges(std::move(payload_ranges));
}

auto DenseReader::validate_payload_ranges(
    std::vector<std::pair<std::uint64_t, std::uint64_t>> ranges) const -> void
{
    std::ranges::sort(ranges);
    for (std::size_t index = 1; index < ranges.size(); ++index)
    {
        if (ranges[index - 1].second > ranges[index].first)
        {
            throw GelexException(
                fmt::format("{}: payloads overlap", path_.string()));
        }
    }
}

auto DenseReader::contains(std::string_view identifier) const -> bool
{
    return index_.contains(identifier);
}

auto DenseReader::size() const noexcept -> std::size_t
{
    return payloads_.size();
}

auto DenseReader::info(std::string_view identifier) const -> const MatrixHeader&
{
    return find_entry(identifier).header;
}

auto DenseReader::payloads() const -> std::vector<MatrixHeader>
{
    std::vector<MatrixHeader> result;
    result.reserve(payloads_.size());
    for (const auto& entry : payloads_)
    {
        result.push_back(entry.header);
    }
    std::ranges::sort(result, {}, &MatrixHeader::identifier);
    return result;
}

auto DenseReader::find_entry(std::string_view identifier) const
    -> const detail::MatrixEntry&
{
    const auto iterator = index_.find(identifier);
    if (iterator == index_.end())
    {
        throw GelexException(
            fmt::format(
                "{}: payload not found: \"{}\"", path_.string(), identifier));
    }
    return payloads_[iterator->second];
}

auto DenseReader::payload_bytes(const detail::MatrixEntry& entry) const
    -> std::span<const std::byte>
{
    const std::span<const std::byte> file{mmap_.data(), mmap_.size()};
    return file.subspan(
        static_cast<std::size_t>(entry.offset),
        static_cast<std::size_t>(entry.size));
}

}  // namespace gelex
