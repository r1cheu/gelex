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

#ifndef GELEX_IO_BINARY_WRITER_H_
#define GELEX_IO_BINARY_WRITER_H_

#include <cassert>
#include <cstddef>
#include <cstdint>
#include <span>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

#include "gelex/exception.h"
#include "gelex/io/binary_format.h"
#include "gelex/io/detail/atomic_output_stream.h"
#include "gelex/io/detail/binary_wire.h"

namespace gelex
{

class BinaryWriter;

template <detail::SupportedDtype T>
class PayloadWriter
{
   public:
    PayloadWriter(const PayloadWriter&) = delete;
    PayloadWriter(PayloadWriter&& other) noexcept;
    auto operator=(const PayloadWriter&) -> PayloadWriter& = delete;
    auto operator=(PayloadWriter&& other) noexcept -> PayloadWriter&;
    ~PayloadWriter() noexcept = default;

    auto append(std::span<const T> payload) -> void;
    auto append(T value) -> void;
    auto write(std::span<const T> payload) -> void;

    [[nodiscard]] auto rows() const noexcept -> std::uint64_t
    {
        return shape_[0];
    }

    [[nodiscard]] auto identifier() const noexcept -> std::string_view
    {
        return identifier_;
    }

   private:
    friend class BinaryWriter;

    PayloadWriter(
        BinaryWriter* writer,
        std::size_t index,
        BinaryShape shape,
        std::string identifier) noexcept
        : writer_(writer),
          index_(index),
          shape_(shape),
          identifier_(std::move(identifier))
    {
    }

    BinaryWriter* writer_;
    std::size_t index_;
    BinaryShape shape_;
    std::string identifier_;
    std::uint64_t columns_written_{};
};

class BinaryWriter
{
   public:
    explicit BinaryWriter(std::string_view output_path);

    BinaryWriter(const BinaryWriter&) = delete;
    BinaryWriter(BinaryWriter&&) = delete;
    auto operator=(const BinaryWriter&) -> BinaryWriter& = delete;
    auto operator=(BinaryWriter&&) -> BinaryWriter& = delete;
    ~BinaryWriter() noexcept;

    template <detail::SupportedDtype T>
    [[nodiscard]] auto reserve(
        std::string_view identifier,
        BinaryShape shape) & -> PayloadWriter<T>
    {
        auto owned_identifier = std::string{identifier};
        const auto index = reserve_payload(
            owned_identifier, detail::binary_type_for<T>, shape);
        return PayloadWriter<T>{
            this, index, shape, std::move(owned_identifier)};
    }

   private:
    template <detail::SupportedDtype>
    friend class PayloadWriter;

    struct PayloadReservation
    {
        detail::PayloadEntry entry;
        std::uint64_t cursor{0};
    };

    auto reserve_payload(
        std::string_view identifier,
        BinaryType type,
        BinaryShape shape) -> std::size_t;
    auto append_bytes(std::size_t index, std::span<const std::byte> bytes)
        -> void;
    auto finalize() -> void;
    auto check_duplicate_identifier(std::string_view identifier) const -> void;
    static auto align_up(std::uint64_t value, std::uint64_t alignment)
        -> std::uint64_t;

    auto write_payload_entry(const detail::PayloadEntry& entry) -> void;
    auto write_footer(
        std::uint64_t directory_offset,
        std::uint64_t payload_count) -> void;

    std::vector<PayloadReservation> reservations_;

    detail::AtomicOutputStream file_;
    std::uint64_t next_offset_{0};
    std::uint64_t file_cursor_{0};
};

template <detail::SupportedDtype T>
PayloadWriter<T>::PayloadWriter(PayloadWriter&& other) noexcept
    : writer_(std::exchange(other.writer_, nullptr)),
      index_(other.index_),
      shape_(other.shape_),
      identifier_(std::move(other.identifier_)),
      columns_written_(other.columns_written_)
{
}

template <detail::SupportedDtype T>
auto PayloadWriter<T>::operator=(PayloadWriter&& other) noexcept
    -> PayloadWriter&
{
    if (this != &other)
    {
        writer_ = std::exchange(other.writer_, nullptr);
        index_ = other.index_;
        shape_ = other.shape_;
        identifier_ = std::move(other.identifier_);
        columns_written_ = other.columns_written_;
    }
    return *this;
}

template <detail::SupportedDtype T>
auto PayloadWriter<T>::append(std::span<const T> payload) -> void
{
    if (writer_ == nullptr)
    {
        throw GelexException("payload writer is invalid");
    }
    assert(std::cmp_equal(payload.size(), shape_[0]));
    writer_->append_bytes(index_, std::as_bytes(payload));
    ++columns_written_;
}

template <detail::SupportedDtype T>
auto PayloadWriter<T>::append(T value) -> void
{
    append(std::span<const T>{&value, 1});
}

template <detail::SupportedDtype T>
auto PayloadWriter<T>::write(std::span<const T> payload) -> void
{
    if (writer_ == nullptr)
    {
        throw GelexException("payload writer is invalid");
    }
    assert(columns_written_ == 0);
    assert(std::cmp_equal(payload.size(), shape_[0] * shape_[1]));
    writer_->append_bytes(index_, std::as_bytes(payload));
    columns_written_ = shape_[1];
}

}  // namespace gelex

#endif  // GELEX_IO_BINARY_WRITER_H_
