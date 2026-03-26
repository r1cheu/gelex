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
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <string_view>
#include <type_traits>
#include <vector>

#include "gelex/io/binary_format.h"

namespace gelex::detail
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

    template <typename T, std::integral Rows, std::integral Cols>
        requires std::is_arithmetic_v<T>
    auto reserve(std::string_view path, Rows rows, Cols cols)
        -> SectionHandle<T>
    {
        return {reserve(
            path,
            binary_format::kTypeByte<T>,
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

    template <MatrixBuffer Matrix>
    auto write(std::string_view path, const Matrix& data) -> void
    {
        using T = element_t<Matrix>;
        auto h = reserve<T>(path, data.rows(), data.cols());
        write(h, data);
    }

    auto finalize() -> void;

   private:
    struct ReservedSection
    {
        TocEntry entry;
        uint64_t cursor{0};
    };

    auto
    reserve(std::string_view path, uint8_t dtype, uint64_t rows, uint64_t cols)
        -> size_t;
    auto write_raw(size_t handle, const char* data, std::streamsize bytes)
        -> void;
    auto check_duplicate_path(std::string_view path) const -> void;
    static auto align_up(uint64_t value, uint64_t alignment) noexcept
        -> uint64_t;
    auto write_footer(uint64_t toc_offset, uint64_t n_sections) -> void;

    std::filesystem::path output_path_;
    std::vector<ReservedSection> reserved_;

    std::ofstream file_;
    uint64_t next_offset_{0};
    uint64_t file_cursor_{0};
};

}  // namespace gelex::detail

#endif  // GELEX_IO_BINARY_WRITER_H_
