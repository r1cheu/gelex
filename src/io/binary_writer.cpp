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

#include <cassert>
#include <concepts>
#include <cstddef>
#include <cstdint>
#include <fmt/format.h>
#include <ios>
#include <limits>
#include <span>
#include <string>
#include <string_view>
#include <utility>

#include "gelex/exception.h"
#include "gelex/infra/log.h"
#include "gelex/io/detail/binary_wire.h"

namespace gelex
{

namespace
{

template <std::unsigned_integral T>
auto write_integer(detail::AtomicOutputStream& file, T value) -> void
{
    file.write(
        reinterpret_cast<const char*>(&value),
        static_cast<std::streamsize>(sizeof(value)));
}

}  // namespace

BinaryWriter::BinaryWriter(std::string_view output_path)
    : file_(std::string(output_path), std::ios::binary | std::ios::trunc)
{
}

BinaryWriter::~BinaryWriter() noexcept
{
    if (closed_ || reservations_.empty())
    {
        return;
    }
    warn(
        fmt::format(
            "{}: destroyed without close(), discarding {} reserved payload(s)",
            file_.path().string(),
            reservations_.size()));
}

auto BinaryWriter::close() -> void
{
    if (closed_)
    {
        throw GelexException(
            fmt::format("{}: writer is closed", file_.path().string()));
    }

    // A reservation is an upper bound; the directory records what was actually
    // written so a shortfall costs the tail, not the whole file.
    std::size_t truncated = 0;
    for (auto& reservation : reservations_)
    {
        auto& entry = reservation.entry;
        const auto written = reservation.cursor - entry.offset;
        if (written == entry.size)
        {
            continue;
        }

        auto& descriptor = entry.info.descriptor;
        const auto column_bytes = detail::checked_product(
            descriptor.shape[0], detail::binary_type_size(descriptor.type));
        if (written % column_bytes != 0)
        {
            throw GelexException(
                fmt::format(
                    "{}: payload \"{}\" ends mid-column: {} of {} bytes "
                    "written",
                    file_.path().string(),
                    entry.info.identifier,
                    written,
                    entry.size));
        }
        descriptor.shape[1] = written / column_bytes;
        entry.size = written;
        ++truncated;
    }
    if (truncated != 0)
    {
        warn(
            fmt::format(
                "{}: {} of {} payload(s) partially written, committing the "
                "columns written so far",
                file_.path().string(),
                truncated,
                reservations_.size()));
    }

    const auto directory_offset
        = align_up(next_offset_, detail::payload_alignment);
    file_.seek(static_cast<std::streamoff>(directory_offset));
    for (const auto& reservation : reservations_)
    {
        write_payload_entry(reservation.entry);
    }
    write_footer(
        directory_offset, static_cast<std::uint64_t>(reservations_.size()));
    file_.commit();
    closed_ = true;
}

auto BinaryWriter::check_duplicate_identifier(std::string_view identifier) const
    -> void
{
    for (const auto& reservation : reservations_)
    {
        if (reservation.entry.info.identifier == identifier)
        {
            throw GelexException(
                fmt::format(
                    "{}: duplicate payload identifier \"{}\"",
                    file_.path().string(),
                    identifier));
        }
    }
}

auto BinaryWriter::reserve_payload(
    std::string_view identifier,
    BinaryType type,
    BinaryShape shape) -> std::size_t
{
    if (closed_)
    {
        throw GelexException(
            fmt::format("{}: writer is closed", file_.path().string()));
    }
    if (identifier.empty())
    {
        throw GelexException(
            fmt::format(
                "{}: payload identifier is empty", file_.path().string()));
    }
    if (identifier.size() > std::numeric_limits<std::uint32_t>::max())
    {
        throw GelexException(
            fmt::format(
                "{}: payload identifier is too long", file_.path().string()));
    }

    check_duplicate_identifier(identifier);

    const auto elements = detail::checked_product(shape[0], shape[1]);
    const auto bytes
        = detail::checked_product(elements, detail::binary_type_size(type));

    const auto aligned_offset
        = align_up(next_offset_, detail::payload_alignment);
    detail::PayloadEntry entry{
        .info
        = PayloadInfo{.identifier = std::string{identifier}, .descriptor = PayloadDescriptor{.type = type, .shape = shape}},
        .offset = aligned_offset,
        .size = bytes};
    next_offset_ = detail::checked_add(aligned_offset, bytes);

    reservations_.push_back(
        PayloadReservation{
            .entry = std::move(entry), .cursor = aligned_offset});

    return reservations_.size() - 1;
}

auto BinaryWriter::append_bytes(
    std::size_t index,
    BinaryType type,
    std::span<const std::byte> bytes) -> void
{
    if (closed_)
    {
        throw GelexException(
            fmt::format("{}: writer is closed", file_.path().string()));
    }
    assert(index < reservations_.size());

    if (!std::in_range<std::streamsize>(bytes.size()))
    {
        throw GelexException(
            fmt::format(
                "{}: write size is out of range", file_.path().string()));
    }

    auto& reservation = reservations_[index];
    assert(reservation.entry.info.descriptor.type == type);
    const auto byte_count = static_cast<std::uint64_t>(bytes.size());
    const auto written = reservation.cursor - reservation.entry.offset;
    if (byte_count > reservation.entry.size - written)
    {
        throw GelexException(
            fmt::format(
                "{}: write overflow: cursor={}, bytes={}, limit={}",
                file_.path().string(),
                reservation.cursor,
                byte_count,
                reservation.entry.offset + reservation.entry.size));
    }

    if (reservation.cursor != file_cursor_)
    {
        file_.seek(static_cast<std::streamoff>(reservation.cursor));
    }
    file_.write(
        reinterpret_cast<const char*>(bytes.data()),
        static_cast<std::streamsize>(bytes.size()));
    reservation.cursor += byte_count;
    file_cursor_ = reservation.cursor;
}

auto BinaryWriter::align_up(std::uint64_t value, std::uint64_t alignment)
    -> std::uint64_t
{
    const auto remainder = value % alignment;
    return remainder == 0 ? value
                          : detail::checked_add(value, alignment - remainder);
}

auto BinaryWriter::write_payload_entry(const detail::PayloadEntry& entry)
    -> void
{
    const auto& info = entry.info;
    write_integer(file_, static_cast<std::uint32_t>(info.identifier.size()));
    file_.write(
        info.identifier.data(),
        static_cast<std::streamsize>(info.identifier.size()));

    const auto encoded_type
        = static_cast<std::byte>(std::to_underlying(info.descriptor.type));
    file_.write(reinterpret_cast<const char*>(&encoded_type), 1);

    write_integer(file_, info.descriptor.shape[0]);
    write_integer(file_, info.descriptor.shape[1]);
    write_integer(file_, entry.offset);
    write_integer(file_, entry.size);
}

auto BinaryWriter::write_footer(
    std::uint64_t directory_offset,
    std::uint64_t payload_count) -> void
{
    file_.write(
        reinterpret_cast<const char*>(detail::binary_format_magic.data()),
        static_cast<std::streamsize>(detail::binary_format_magic.size()));
    write_integer(file_, directory_offset);
    write_integer(file_, payload_count);
}

}  // namespace gelex
