// Copyright 2026 RuLei Chen
// SPDX-License-Identifier: Apache-2.0

#ifndef GELEX_IO_CSC_WRITER_H_
#define GELEX_IO_CSC_WRITER_H_

#include <cstddef>
#include <cstdint>
#include <fmt/format.h>
#include <ranges>
#include <span>
#include <string_view>
#include <utility>
#include <vector>

#include "gelex/exception.h"
#include "gelex/io/binary_format.h"
#include "gelex/io/detail/atomic_output_stream.h"

namespace gelex
{

class CscWriter;

// Converts dense columns to CSC entries and hands them to the owning writer,
// which tracks column counts, spools and completion.
template <detail::SupportedDtype T>
class CscStream
{
   public:
    CscStream(const CscStream&) = delete;
    CscStream(CscStream&& other) noexcept
        : writer_(std::exchange(other.writer_, nullptr)),
          index_(other.index_),
          shape_(other.shape_),
          values_(std::move(other.values_)),
          indices_(std::move(other.indices_))
    {
    }
    auto operator=(const CscStream&) -> CscStream& = delete;
    auto operator=(CscStream&& other) noexcept -> CscStream&
    {
        if (this != &other)
        {
            writer_ = std::exchange(other.writer_, nullptr);
            index_ = other.index_;
            shape_ = other.shape_;
            values_ = std::move(other.values_);
            indices_ = std::move(other.indices_);
        }
        return *this;
    }
    ~CscStream() noexcept = default;

    // Writes one complete column, omitting exact zeros without scalar
    // conversion.
    auto operator<<(std::span<const T> column) -> CscStream&;

   private:
    friend class CscWriter;

    CscStream(CscWriter* writer, std::size_t index, BinaryShape shape) noexcept
        : writer_(writer), index_(index), shape_(shape)
    {
    }

    CscWriter* writer_;
    std::size_t index_;
    BinaryShape shape_;
    std::vector<T> values_;
    std::vector<std::int64_t> indices_;
};

// Streams borrow this object's address and must not outlive it.
class CscWriter
{
   public:
    explicit CscWriter(std::string_view output_path);

    CscWriter(const CscWriter&) = delete;
    CscWriter(CscWriter&&) = delete;
    auto operator=(const CscWriter&) -> CscWriter& = delete;
    auto operator=(CscWriter&&) -> CscWriter& = delete;
    ~CscWriter() noexcept;

    template <detail::SupportedDtype T>
    [[nodiscard]] auto reserve(
        std::string_view identifier,
        BinaryShape shape) & -> CscStream<T>
    {
        const auto index
            = reserve(identifier, detail::binary_type_for<T>, shape);
        return CscStream<T>{this, index, shape};
    }

    // Requires every matrix to be complete. Only close() publishes the file;
    // destruction discards uncommitted output, including after a failed close,
    // and logs an error unless an exception is unwinding.
    auto close() -> void;
    auto is_open() const noexcept -> bool { return !closed_; }

   private:
    template <detail::SupportedDtype>
    friend class CscStream;

    struct MatrixStorage;

    auto reserve(
        std::string_view identifier,
        BinaryType type,
        BinaryShape shape) -> std::size_t;
    auto append_column(
        std::size_t index,
        std::span<const std::byte> values,
        std::span<const std::int64_t> indices) -> void;

    auto throw_if_closed() const -> void;
    auto validate_identifier(std::string_view identifier) const -> void;
    auto validate_complete() const -> void;
    auto complete(std::size_t index) -> void;
    auto merge(std::size_t index) -> void;
    auto discard_spools() noexcept -> void;
    auto write_index(std::uint64_t offset) -> void;
    auto write_index_entry(const MatrixStorage& matrix) -> void;
    auto write_footer(std::uint64_t index_offset) -> void;

    detail::AtomicOutputStream output_;
    std::vector<MatrixStorage> matrices_;
    std::vector<char> copy_buffer_;
    std::uint64_t offset_{};
    bool closed_{false};
};

[[nodiscard]] auto open_csc_writer(std::string_view output_path) -> CscWriter;

template <detail::SupportedDtype T>
auto CscStream<T>::operator<<(std::span<const T> column) -> CscStream&
{
    if (writer_ == nullptr)
    {
        throw GelexException("CSC stream is invalid");
    }
    if (!std::cmp_equal(column.size(), shape_[0]))
    {
        throw GelexException(
            fmt::format(
                "CSC matrix {}: expected {} rows, got {}",
                index_,
                shape_[0],
                column.size()));
    }
    values_.clear();
    indices_.clear();
    for (auto [row, value] : std::views::enumerate(column))
    {
        if (value != T{0})
        {
            values_.push_back(value);
            indices_.push_back(static_cast<std::int64_t>(row));
        }
    }
    writer_->append_column(index_, std::as_bytes(std::span{values_}), indices_);
    return *this;
}

}  // namespace gelex

#endif  // GELEX_IO_CSC_WRITER_H_
