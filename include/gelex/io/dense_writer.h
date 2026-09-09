// Copyright 2026 RuLei Chen
// SPDX-License-Identifier: Apache-2.0

#ifndef GELEX_IO_DENSE_WRITER_H_
#define GELEX_IO_DENSE_WRITER_H_

#include <cstddef>
#include <cstdint>
#include <fmt/format.h>
#include <span>
#include <string_view>
#include <utility>
#include <vector>

#include "gelex/exception.h"
#include "gelex/io/binary_format.h"
#include "gelex/io/detail/atomic_output_stream.h"

namespace gelex
{

class DenseWriter;

// Hands whole columns of one reserved matrix to the owning writer, which
// tracks the write cursor and the matrix identity.
template <detail::SupportedDtype T>
class DenseStream
{
   public:
    DenseStream(const DenseStream&) = delete;
    DenseStream(DenseStream&& other) noexcept
        : writer_(std::exchange(other.writer_, nullptr)),
          index_(other.index_),
          shape_(other.shape_)
    {
    }
    auto operator=(const DenseStream&) -> DenseStream& = delete;
    auto operator=(DenseStream&& other) noexcept -> DenseStream&
    {
        if (this != &other)
        {
            writer_ = std::exchange(other.writer_, nullptr);
            index_ = other.index_;
            shape_ = other.shape_;
        }
        return *this;
    }
    ~DenseStream() noexcept = default;

    // Appends one complete column or the whole column-major matrix.
    auto operator<<(std::span<const T> values) -> DenseStream&;
    auto operator<<(T value) -> DenseStream&
    {
        return *this << std::span<const T>{&value, 1};
    }

    [[nodiscard]] auto identifier() const -> std::string_view;

   private:
    friend class DenseWriter;

    DenseStream(
        DenseWriter* writer,
        std::size_t index,
        BinaryShape shape) noexcept
        : writer_(writer), index_(index), shape_(shape)
    {
    }

    auto throw_if_invalid() const -> void;

    DenseWriter* writer_;
    std::size_t index_;
    BinaryShape shape_;
};

// Streams borrow this object's address and must not outlive it.
class DenseWriter
{
   public:
    explicit DenseWriter(std::string_view output_path);

    DenseWriter(const DenseWriter&) = delete;
    DenseWriter(DenseWriter&&) = delete;
    auto operator=(const DenseWriter&) -> DenseWriter& = delete;
    auto operator=(DenseWriter&&) -> DenseWriter& = delete;
    ~DenseWriter() noexcept;

    template <detail::SupportedDtype T>
    [[nodiscard]] auto reserve(
        std::string_view identifier,
        BinaryShape shape) & -> DenseStream<T>
    {
        const auto index
            = reserve(identifier, detail::binary_type_for<T>, shape);
        return DenseStream<T>{this, index, shape};
    }

    // Requires every matrix to be complete. Only close() publishes the
    // file; destruction discards uncommitted output, including after a
    // failed close, and logs an error unless an exception is unwinding.
    auto close() -> void;
    [[nodiscard]] auto is_open() const noexcept -> bool { return !closed_; }

   private:
    template <detail::SupportedDtype>
    friend class DenseStream;

    struct Reservation;

    auto reserve(
        std::string_view identifier,
        BinaryType type,
        BinaryShape shape) -> std::size_t;
    auto append_bytes(std::size_t index, std::span<const std::byte> bytes)
        -> void;
    [[nodiscard]] auto identifier(std::size_t index) const -> std::string_view;

    auto throw_if_closed() const -> void;
    auto validate_identifier(std::string_view identifier) const -> void;
    auto validate_complete() const -> void;
    auto write_index(std::uint64_t offset) -> void;
    auto write_footer(std::uint64_t index_offset) -> void;

    detail::AtomicOutputStream output_;
    std::vector<Reservation> reservations_;
    std::uint64_t next_offset_{};
    bool closed_{false};
};

[[nodiscard]] auto open_dense_writer(std::string_view output_path)
    -> DenseWriter;

template <detail::SupportedDtype T>
auto DenseStream<T>::throw_if_invalid() const -> void
{
    if (writer_ == nullptr)
    {
        throw GelexException("dense stream is invalid");
    }
}

template <detail::SupportedDtype T>
auto DenseStream<T>::identifier() const -> std::string_view
{
    throw_if_invalid();
    return writer_->identifier(index_);
}

template <detail::SupportedDtype T>
auto DenseStream<T>::operator<<(std::span<const T> values) -> DenseStream&
{
    throw_if_invalid();
    if (!std::cmp_equal(values.size(), shape_[0])
        && !std::cmp_equal(values.size(), shape_[0] * shape_[1]))
    {
        throw GelexException(
            fmt::format(
                "dense matrix {}: expected {} column values or {} matrix "
                "values, got {}",
                index_,
                shape_[0],
                shape_[0] * shape_[1],
                values.size()));
    }
    writer_->append_bytes(index_, std::as_bytes(values));
    return *this;
}

}  // namespace gelex

#endif  // GELEX_IO_DENSE_WRITER_H_
