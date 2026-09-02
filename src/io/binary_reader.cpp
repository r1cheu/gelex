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

#include "gelex/io/binary_reader.h"

#include <algorithm>
#include <array>
#include <cstddef>
#include <cstdint>
#include <cstring>
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

class BinaryCursor
{
   public:
    BinaryCursor(
        std::span<const std::byte> bytes,
        std::string_view path) noexcept
        : remaining_(bytes), path_(path)
    {
    }

    auto read_bytes(std::size_t size, std::string_view field)
        -> std::span<const std::byte>
    {
        if (size > remaining_.size())
        {
            throw GelexException(fmt::format("{}: truncated {}", path_, field));
        }
        const auto bytes = remaining_.first(size);
        remaining_ = remaining_.subspan(size);
        return bytes;
    }

    [[nodiscard]] auto size() const noexcept -> std::size_t
    {
        return remaining_.size();
    }

    [[nodiscard]] auto empty() const noexcept -> bool
    {
        return remaining_.empty();
    }

   private:
    std::span<const std::byte> remaining_;
    std::string_view path_;
};

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

    const auto footer = file.last(detail::footer_size);
    if (!std::ranges::equal(
            detail::binary_format_magic,
            footer.first(detail::binary_format_magic.size())))
    {
        throw GelexException(fmt::format("{}: invalid container magic", path));
    }

    std::array<std::uint64_t, 2> footer_values{};
    std::memcpy(
        footer_values.data(),
        footer.subspan(detail::binary_format_magic.size()).data(),
        sizeof(footer_values));
    const auto directory_offset = footer_values[0];
    const auto payload_count = footer_values[1];
    const auto footer_offset = file.size() - detail::footer_size;
    if (directory_offset % detail::payload_alignment != 0
        || directory_offset > footer_offset)
    {
        throw GelexException(fmt::format("{}: invalid directory offset", path));
    }

    return Footer{
        .directory_offset = directory_offset, .payload_count = payload_count};
}

auto parse_payload_entry(
    BinaryCursor& cursor,
    const Footer& footer,
    std::uint64_t index,
    std::string_view path) -> detail::PayloadEntry
{
    const auto identifier_size_bytes
        = cursor.read_bytes(sizeof(std::uint32_t), "payload entry");
    const auto identifier_size
        = detail::read_integer<std::uint32_t>(identifier_size_bytes.data());
    const auto identifier_bytes
        = cursor.read_bytes(identifier_size, "payload identifier");
    std::string identifier{
        reinterpret_cast<const char*>(identifier_bytes.data()),
        identifier_bytes.size()};

    const auto encoded_type
        = cursor.read_bytes(1, "payload descriptor").front();
    const auto type = detail::decode_binary_type(encoded_type);
    const auto descriptor_bytes
        = cursor.read_bytes(4 * sizeof(std::uint64_t), "payload descriptor");
    std::array<std::uint64_t, 4> descriptor_values{};
    std::memcpy(
        descriptor_values.data(),
        descriptor_bytes.data(),
        sizeof(descriptor_values));
    const BinaryShape shape{descriptor_values[0], descriptor_values[1]};
    const auto offset = descriptor_values[2];
    const auto size = descriptor_values[3];

    if (offset % detail::payload_alignment != 0
        || offset > footer.directory_offset
        || size > footer.directory_offset - offset)
    {
        throw GelexException(
            fmt::format(
                "{}: payload {} is outside the data region", path, index));
    }

    const auto element_count = detail::checked_product(shape[0], shape[1]);
    if (size
        != detail::checked_product(
            element_count, detail::binary_type_size(type)))
    {
        throw GelexException(
            fmt::format(
                "{}: payload {} size does not match shape", path, index));
    }

    PayloadInfo info{
        .identifier = std::move(identifier),
        .descriptor = PayloadDescriptor{.type = type, .shape = shape}};
    return detail::PayloadEntry{
        .info = std::move(info), .offset = offset, .size = size};
}

}  // namespace

BinaryReader::BinaryReader(std::string_view file_path)
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

    parse_footer_and_directory();
}

auto BinaryReader::parse_footer_and_directory() -> void
{
    const std::span<const std::byte> file{mmap_.data(), mmap_.size()};
    const auto path_string = path_.string();
    const auto footer = parse_footer(file, path_string);
    const auto footer_offset = file.size() - detail::footer_size;
    const auto directory_offset
        = static_cast<std::size_t>(footer.directory_offset);
    BinaryCursor directory{
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
            = index_.emplace(payload.info.identifier, index_.size());
        if (!inserted)
        {
            throw GelexException(
                fmt::format(
                    "{}: duplicate payload identifier \"{}\"",
                    path_string,
                    payload.info.identifier));
        }
    }

    validate_payload_ranges(std::move(payload_ranges));
}

auto BinaryReader::validate_payload_ranges(
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

auto BinaryReader::contains(std::string_view identifier) const -> bool
{
    return index_.contains(identifier);
}

auto BinaryReader::size() const noexcept -> std::size_t
{
    return payloads_.size();
}

auto BinaryReader::info(std::string_view identifier) const& -> const
    PayloadInfo&
{
    return find_entry(identifier).info;
}

auto BinaryReader::payloads() const -> std::vector<PayloadInfo>
{
    std::vector<PayloadInfo> result;
    result.reserve(payloads_.size());
    for (const auto& entry : payloads_)
    {
        result.push_back(entry.info);
    }
    std::ranges::sort(result, {}, &PayloadInfo::identifier);
    return result;
}

auto BinaryReader::find_entry(std::string_view identifier) const
    -> const detail::PayloadEntry&
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

auto BinaryReader::payload_bytes(const detail::PayloadEntry& entry) const
    -> std::span<const std::byte>
{
    const std::span<const std::byte> file{mmap_.data(), mmap_.size()};
    return file.subspan(
        static_cast<std::size_t>(entry.offset),
        static_cast<std::size_t>(entry.size));
}

}  // namespace gelex
