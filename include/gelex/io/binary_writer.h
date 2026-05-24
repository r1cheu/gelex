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

#include <concepts>
#include <cstddef>
#include <cstdint>
#include <ios>
#include <ranges>
#include <string_view>
#include <type_traits>
#include <utility>
#include <vector>

#include "gelex/io/detail/atomic_ofstream.h"
#include "gelex/io/detail/binary_format.h"

namespace gelex::io
{

template <typename C>
using element_t = std::decay_t<decltype(*std::declval<const C&>().data())>;

template <typename C, typename T = element_t<C>>
concept DataBuffer = std::is_arithmetic_v<T> && requires(const C& c) {
    { c.data() } -> std::convertible_to<const T*>;
    { c.size() } -> std::convertible_to<size_t>;
};

template <typename C, typename T = element_t<C>>
concept MatrixBuffer = DataBuffer<C, T> && requires(const C& c) {
    { c.rows() } -> std::convertible_to<size_t>;
    { c.cols() } -> std::convertible_to<size_t>;
};

template <typename T>
struct SectionHandle
{
    size_t index;
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

    template <typename T, std::integral Rows, std::integral Cols>
        requires std::is_arithmetic_v<T>
    auto reserve(std::string_view path, Rows rows, Cols cols)
        -> SectionHandle<T>
    {
        return {reserve(
            path,
            detail::kTypeByte<T>,
            static_cast<uint64_t>(rows),
            static_cast<uint64_t>(cols))};
    }

    template <typename T, DataBuffer<T> Container>
    auto write(SectionHandle<T> handle, const Container& data) -> void
    {
        write_raw(
            handle.index,
            reinterpret_cast<const char*>(data.data()),
            static_cast<std::streamsize>(
                static_cast<size_t>(data.size()) * sizeof(T)));
    }

    template <typename T>
        requires std::is_arithmetic_v<T>
    auto write(SectionHandle<T> handle, T value) -> void
    {
        write_raw(
            handle.index,
            reinterpret_cast<const char*>(&value),
            static_cast<std::streamsize>(sizeof(T)));
    }

    template <typename T>
        requires std::is_arithmetic_v<T>
    auto write(std::string_view path, T value) -> void
    {
        auto h = reserve<T>(path, 1, 1);
        write(h, value);
    }

    template <MatrixBuffer Matrix>
    auto write(std::string_view path, const Matrix& data) -> void
    {
        using T = element_t<Matrix>;
        auto h = reserve<T>(path, data.rows(), data.cols());
        write(h, data);
    }

    template <std::ranges::forward_range R>
        requires std::
            convertible_to<std::ranges::range_value_t<R>, std::string_view>
        auto write_strings(std::string_view path, const R& names) -> void
    {
        uint64_t total_bytes = 0;
        for (const std::string_view s : names)
        {
            total_bytes += s.size() + 1;
        }

        const auto handle = reserve_strings(
            path,
            static_cast<uint64_t>(std::ranges::distance(names)),
            total_bytes);

        for (const std::string_view s : names)
        {
            write_raw(handle, s.data(), static_cast<std::streamsize>(s.size()));
            write_raw(handle, "\0", 1);
        }
    }

   private:
    struct ReservedSection
    {
        detail::TocEntry entry;
        uint64_t cursor{0};
    };

    auto finalize() -> void;
    auto
    reserve(std::string_view path, uint8_t dtype, uint64_t rows, uint64_t cols)
        -> size_t;
    auto reserve_strings(std::string_view path, uint64_t rows, uint64_t bytes)
        -> size_t;
    auto reserve_section(
        std::string_view path,
        uint8_t dtype,
        uint64_t rows,
        uint64_t cols,
        uint64_t bytes) -> size_t;
    auto write_raw(size_t handle, const char* data, std::streamsize bytes)
        -> void;
    auto check_duplicate_path(std::string_view path) const -> void;
    static auto align_up(uint64_t value, uint64_t alignment) noexcept
        -> uint64_t;
    auto write_footer(uint64_t toc_offset, uint64_t n_sections) -> void;

    std::vector<ReservedSection> reserved_;

    detail::AtomicOfstream file_;
    uint64_t next_offset_{0};
    uint64_t file_cursor_{0};
};

}  // namespace gelex::io

#endif  // GELEX_IO_BINARY_WRITER_H_
