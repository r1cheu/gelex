// Copyright 2026 RuLei Chen
// SPDX-License-Identifier: Apache-2.0

#include "gelex/io/dense_writer.h"

#include <cstddef>
#include <cstdint>
#include <exception>
#include <fmt/format.h>
#include <ios>
#include <span>
#include <string>
#include <string_view>
#include <utility>

#include "gelex/exception.h"
#include "gelex/infra/log.h"
#include "gelex/io/binary_format.h"
#include "gelex/io/detail/atomic_output_stream.h"
#include "gelex/io/detail/binary_wire.h"

namespace gelex
{

struct DenseWriter::Reservation
{
    MatrixHeader header;
    std::uint64_t offset{};
    std::uint64_t size{};
    std::uint64_t cursor{};
};

DenseWriter::DenseWriter(std::string_view output_path)
    : output_(std::string{output_path})
{
}

DenseWriter::~DenseWriter() noexcept
{
    if (!closed_ && std::uncaught_exceptions() == 0)
    {
        try
        {
            error(
                fmt::format(
                    "{}: unclosed dense writer destroyed, output "
                    "discarded",
                    output_.path().string()));
        }
        catch (...)  // NOLINT(bugprone-empty-catch): dtor must be noexcept
        {
        }
    }
}

auto open_dense_writer(std::string_view output_path) -> DenseWriter
{
    return DenseWriter{output_path};
}

auto DenseWriter::throw_if_closed() const -> void
{
    if (closed_)
    {
        throw GelexException(
            fmt::format("{}: dense writer is closed", output_.path().string()));
    }
}

auto DenseWriter::validate_identifier(std::string_view identifier) const -> void
{
    if (identifier.empty())
    {
        throw GelexException(
            fmt::format(
                "{}: matrix identifier is empty", output_.path().string()));
    }
    for (const auto& reservation : reservations_)
    {
        if (reservation.header.identifier == identifier)
        {
            throw GelexException(
                fmt::format(
                    "{}: duplicate matrix identifier \"{}\"",
                    output_.path().string(),
                    identifier));
        }
    }
}

auto DenseWriter::identifier(std::size_t index) const -> std::string_view
{
    return reservations_[index].header.identifier;
}

auto DenseWriter::reserve(
    std::string_view identifier,
    BinaryType type,
    BinaryShape shape) -> std::size_t
{
    throw_if_closed();
    validate_identifier(identifier);
    const auto bytes = shape[0] * shape[1] * detail::binary_type_size(type);
    const auto offset = detail::align_payload_offset(next_offset_);
    const auto end = offset + bytes;
    const auto index = reservations_.size();
    reservations_.push_back(
        Reservation{
            .header
            = {.identifier = std::string{identifier},
               .type = type,
               .shape = shape},
            .offset = offset,
            .size = bytes,
            .cursor = offset});
    next_offset_ = end;
    return index;
}

auto DenseWriter::append_bytes(
    std::size_t index,
    std::span<const std::byte> bytes) -> void
{
    throw_if_closed();
    auto& reservation = reservations_[index];
    const auto written = reservation.cursor - reservation.offset;
    if (bytes.size() > reservation.size - written)
    {
        throw GelexException(
            fmt::format(
                "{}: matrix \"{}\" overflow: {} of {} bytes written, {} "
                "more requested",
                output_.path().string(),
                reservation.header.identifier,
                written,
                reservation.size,
                bytes.size()));
    }
    output_.seek(static_cast<std::streamoff>(reservation.cursor));
    output_.write(
        reinterpret_cast<const char*>(bytes.data()),
        static_cast<std::streamsize>(bytes.size()));
    reservation.cursor += bytes.size();
}

auto DenseWriter::validate_complete() const -> void
{
    for (const auto& reservation : reservations_)
    {
        const auto written = reservation.cursor - reservation.offset;
        if (written != reservation.size)
        {
            throw GelexException(
                fmt::format(
                    "{}: matrix \"{}\" is incomplete: {} of {} bytes written",
                    output_.path().string(),
                    reservation.header.identifier,
                    written,
                    reservation.size));
        }
    }
}

auto DenseWriter::write_index(std::uint64_t offset) -> void
{
    auto& output = output_;
    output.seek(static_cast<std::streamoff>(offset));
    for (const auto& reservation : reservations_)
    {
        const auto& header = reservation.header;
        detail::write_integer(
            output, static_cast<std::uint32_t>(header.identifier.size()));
        output.write(header.identifier);
        detail::write_integer(output, std::to_underlying(header.type));
        detail::write_integer(output, header.shape[0]);
        detail::write_integer(output, header.shape[1]);
        detail::write_integer(output, reservation.offset);
        detail::write_integer(output, reservation.size);
    }
}

auto DenseWriter::write_footer(std::uint64_t index_offset) -> void
{
    auto& output = output_;
    output.write(
        reinterpret_cast<const char*>(detail::binary_format_magic.data()),
        static_cast<std::streamsize>(detail::binary_format_magic.size()));
    detail::write_integer(output, index_offset);
    detail::write_integer(
        output, static_cast<std::uint64_t>(reservations_.size()));
}

auto DenseWriter::close() -> void
{
    if (closed_)
    {
        return;
    }
    validate_complete();
    const auto index_offset = detail::align_payload_offset(next_offset_);
    write_index(index_offset);
    write_footer(index_offset);
    output_.commit();
    closed_ = true;
}

}  // namespace gelex
