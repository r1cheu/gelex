// Copyright 2026 RuLei Chen
// SPDX-License-Identifier: Apache-2.0

#include "gelex/io/csc_writer.h"

#include <algorithm>
#include <array>
#include <cstddef>
#include <cstdint>
#include <exception>
#include <filesystem>
#include <fmt/format.h>
#include <fstream>
#include <ios>
#include <span>
#include <string>
#include <string_view>
#include <system_error>
#include <utility>
#include <vector>

#include "gelex/exception.h"
#include "gelex/infra/log.h"
#include "gelex/io/binary_format.h"
#include "gelex/io/detail/atomic_output_stream.h"
#include "gelex/io/detail/binary_wire.h"

namespace gelex
{

namespace
{

// CSC v1 stores aligned values, int64 indices and int64 indptr arrays.
// Blocks follow completion order; index entries follow reservation order.
// Index entries: name length/name, dtype, rows/cols/nnz, three offsets.
// Footer: magic, uint64 index offset, uint64 matrix count (little endian).
enum class CscArray : std::uint8_t
{
    values,
    indices,
    indptr,
};
constexpr std::array csc_arrays{
    CscArray::values,
    CscArray::indices,
    CscArray::indptr};
constexpr std::size_t array_count = detail::csc_array_count;
static_assert(csc_arrays.size() == array_count);
constexpr std::size_t copy_buffer_bytes = 64ULL * 1024;

auto merge_csc_array(
    const std::filesystem::path& path,
    std::uint64_t bytes,
    detail::AtomicOutputStream& output,
    std::span<char> buffer) -> void
{
    std::ifstream input(path, std::ios::binary);
    if (!input)
    {
        throw GelexException(
            fmt::format("{}: failed to open CSC array", path.string()));
    }
    while (bytes != 0)
    {
        const auto count = static_cast<std::streamsize>(
            std::min(bytes, static_cast<std::uint64_t>(buffer.size())));
        input.read(buffer.data(), count);
        if (!input)
        {
            throw GelexException(
                fmt::format("{}: failed to read CSC array", path.string()));
        }
        output.write(buffer.data(), count);
        bytes -= static_cast<std::uint64_t>(count);
    }
}

// Temporary per-matrix values/indices/indptr files "<output>.N". Each array
// is written through its own AtomicOutputStream, so "<output>.N.tmp" is
// created exclusively and a failed write discards the array; close() commits
// the arrays so merge() can read them back.
class Spool
{
   public:
    Spool(const std::filesystem::path& output, std::size_t matrix_index)
        : streams_{
              open_array(output, matrix_index, CscArray::values),
              open_array(output, matrix_index, CscArray::indices),
              open_array(output, matrix_index, CscArray::indptr)}
    {
    }

    auto write(CscArray array, std::span<const std::byte> bytes) -> void
    {
        stream_of(array).write(
            reinterpret_cast<const char*>(bytes.data()),
            static_cast<std::streamsize>(bytes.size()));
    }

    auto close() -> void
    {
        for (auto& stream : streams_)
        {
            stream.commit();
        }
    }

    // Valid after close().
    [[nodiscard]] auto path(CscArray array) const
        -> const std::filesystem::path&
    {
        return streams_[std::to_underlying(array)].path();
    }

    // Uncommitted arrays are discarded by their streams. Any committed path
    // is ours, because open_array() refused pre-existing ones.
    auto remove() noexcept -> void
    {
        for (const auto& stream : streams_)
        {
            std::error_code ec;
            std::filesystem::remove(stream.path(), ec);
        }
    }

   private:
    static auto open_array(
        const std::filesystem::path& output,
        std::size_t matrix_index,
        CscArray array) -> detail::AtomicOutputStream
    {
        auto path = output;
        path += fmt::format(
            ".{}",
            (matrix_index * array_count) + std::to_underlying(array) + 1);
        // commit() would replace this path; refuse to clobber a file that
        // is not ours.
        if (std::filesystem::exists(path))
        {
            throw GelexException(
                fmt::format(
                    "{}: CSC temporary file already exists", path.string()));
        }
        return detail::AtomicOutputStream{std::move(path)};
    }

    auto stream_of(CscArray array) -> detail::AtomicOutputStream&
    {
        return streams_[std::to_underlying(array)];
    }

    std::array<detail::AtomicOutputStream, array_count> streams_;
};

}  // namespace

struct CscWriter::MatrixStorage
{
    MatrixStorage(MatrixHeader header, Spool spool)
        : header(std::move(header)), spool(std::move(spool))
    {
    }

    MatrixHeader header;
    Spool spool;
    std::array<std::uint64_t, array_count> offsets{};
    std::uint64_t columns_written{};
    std::uint64_t nnz{};
    bool complete{false};
};

CscWriter::CscWriter(std::string_view output_path)
    : output_(std::string{output_path}), copy_buffer_(copy_buffer_bytes)
{
}

CscWriter::~CscWriter() noexcept
{
    discard_spools();
    if (!closed_ && std::uncaught_exceptions() == 0)
    {
        try
        {
            error(
                fmt::format(
                    "{}: unclosed CSC writer destroyed, output "
                    "discarded",
                    output_.path().string()));
        }
        catch (...)  // NOLINT(bugprone-empty-catch): dtor must be noexcept
        {
        }
    }
}

auto CscWriter::discard_spools() noexcept -> void
{
    for (auto& matrix : matrices_)
    {
        matrix.spool.remove();
    }
}

auto open_csc_writer(std::string_view output_path) -> CscWriter
{
    return CscWriter{output_path};
}

auto CscWriter::throw_if_closed() const -> void
{
    if (closed_)
    {
        throw GelexException(
            fmt::format("{}: CSC writer is closed", output_.path().string()));
    }
}

auto CscWriter::validate_identifier(std::string_view identifier) const -> void
{
    if (identifier.empty())
    {
        throw GelexException("CSC identifier is empty");
    }
    for (const auto& matrix : matrices_)
    {
        if (matrix.header.identifier == identifier)
        {
            throw GelexException(
                fmt::format("duplicate CSC identifier \"{}\"", identifier));
        }
    }
}

auto CscWriter::reserve(
    std::string_view identifier,
    BinaryType type,
    BinaryShape shape) -> std::size_t
{
    throw_if_closed();
    validate_identifier(identifier);
    const auto index = matrices_.size();
    MatrixStorage matrix{
        MatrixHeader{
            .identifier = std::string{identifier},
            .type = type,
            .shape = shape},
        Spool{output_.path(), index}};
    const std::int64_t pointer = 0;
    matrix.spool.write(CscArray::indptr, std::as_bytes(std::span{&pointer, 1}));
    matrices_.push_back(std::move(matrix));
    if (shape[1] == 0)
    {
        complete(index);
    }
    return index;
}

auto CscWriter::append_column(
    std::size_t index,
    std::span<const std::byte> values,
    std::span<const std::int64_t> indices) -> void
{
    throw_if_closed();
    auto& matrix = matrices_[index];
    if (matrix.columns_written == matrix.header.shape[1])
    {
        throw GelexException(
            fmt::format(
                "{}: CSC column count exceeded", matrix.header.identifier));
    }
    matrix.spool.write(CscArray::values, values);
    matrix.spool.write(CscArray::indices, std::as_bytes(indices));
    matrix.nnz += indices.size();
    const auto pointer = static_cast<std::int64_t>(matrix.nnz);
    matrix.spool.write(CscArray::indptr, std::as_bytes(std::span{&pointer, 1}));
    if (++matrix.columns_written == matrix.header.shape[1])
    {
        complete(index);
    }
}

auto CscWriter::complete(std::size_t index) -> void
{
    auto& matrix = matrices_[index];
    matrix.spool.close();
    merge(index);
    matrix.spool.remove();
    matrix.complete = true;
}

auto CscWriter::merge(std::size_t index) -> void
{
    auto& matrix = matrices_[index];
    const auto& header = matrix.header;
    const auto sizes
        = detail::csc_array_sizes(header.type, header.shape[1], matrix.nnz);
    auto offset = offset_;
    for (const auto array : csc_arrays)
    {
        const auto size = sizes[std::to_underlying(array)];
        const auto aligned_offset = detail::align_payload_offset(offset);
        offset = aligned_offset + size;
        matrix.offsets[std::to_underlying(array)] = aligned_offset;
        output_.seek(static_cast<std::streamoff>(aligned_offset));
        merge_csc_array(matrix.spool.path(array), size, output_, copy_buffer_);
    }
    offset_ = offset;
}

auto CscWriter::validate_complete() const -> void
{
    for (const auto& matrix : matrices_)
    {
        if (!matrix.complete)
        {
            throw GelexException(
                fmt::format(
                    "{}: CSC matrix is incomplete ({} columns required)",
                    matrix.header.identifier,
                    matrix.header.shape[1]));
        }
    }
}

auto CscWriter::write_index_entry(const MatrixStorage& matrix) -> void
{
    auto& output = output_;
    const auto& header = matrix.header;
    detail::write_integer(
        output, static_cast<std::uint32_t>(header.identifier.size()));
    output.write(header.identifier);
    detail::write_integer(output, std::to_underlying(header.type));
    detail::write_integer(output, header.shape[0]);
    detail::write_integer(output, header.shape[1]);
    detail::write_integer(output, matrix.nnz);
    for (const auto offset : matrix.offsets)
    {
        detail::write_integer(output, offset);
    }
}

auto CscWriter::write_index(std::uint64_t offset) -> void
{
    output_.seek(static_cast<std::streamoff>(offset));
    for (const auto& matrix : matrices_)
    {
        write_index_entry(matrix);
    }
}

auto CscWriter::write_footer(std::uint64_t index_offset) -> void
{
    output_.write(
        reinterpret_cast<const char*>(detail::csc_format_magic.data()),
        static_cast<std::streamsize>(detail::csc_format_magic.size()));
    detail::write_integer(output_, index_offset);
    detail::write_integer(
        output_, static_cast<std::uint64_t>(matrices_.size()));
}

auto CscWriter::close() -> void
{
    if (closed_)
    {
        return;
    }
    validate_complete();
    const auto index_offset = detail::align_payload_offset(offset_);
    write_index(index_offset);
    write_footer(index_offset);
    output_.commit();
    closed_ = true;
}

}  // namespace gelex
