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

#ifndef GELEX_DATA_BINARY_WRITER_H_
#define GELEX_DATA_BINARY_WRITER_H_

#include <array>
#include <cstddef>
#include <cstdint>
#include <filesystem>
#include <format>
#include <fstream>
#include <span>
#include <string>
#include <string_view>
#include <type_traits>
#include <vector>

#include <Eigen/Core>

#include "gelex/exception.h"

#include "parser.h"

namespace gelex::detail
{

template <typename eT>
class BinaryWriter
{
   public:
    static constexpr size_t kDefaultBufferSize = static_cast<size_t>(64 * 1024);
    static constexpr std::array<std::byte, 8> kMagic
        = {std::byte{'G'},
           std::byte{'E'},
           std::byte{'L'},
           std::byte{'E'},
           std::byte{'X'},
           std::byte{'B'},
           std::byte{'W'},
           std::byte{'1'}};
    static constexpr uint32_t kVersion = 1;
    static constexpr size_t kMetaSize = 8 + 4 + 8 + 8 + 1;

    explicit BinaryWriter(std::string_view file_path);

    BinaryWriter(const BinaryWriter&) = delete;
    BinaryWriter(BinaryWriter&&) noexcept = default;
    auto operator=(const BinaryWriter&) -> BinaryWriter& = delete;
    auto operator=(BinaryWriter&&) noexcept -> BinaryWriter& = default;

    ~BinaryWriter() noexcept;

    auto write(
        const Eigen::Ref<const Eigen::Matrix<eT, Eigen::Dynamic, 1>>& record)
        -> void;

    auto finish() -> void;

   private:
    static constexpr auto is_supported_type() -> bool;
    static constexpr auto dtype_code() -> uint8_t;

    auto write_meta() -> void;
    auto write_u32_le(uint32_t value) -> void;
    auto write_u64_le(uint64_t value) -> void;

    std::filesystem::path path_;
    std::vector<char> io_buffer_;
    std::ofstream file_;
    uint64_t n_rows_ = 0;
    uint64_t n_cols_ = 0;
    bool has_shape_ = false;
    bool finished_ = false;
};

template <typename eT>
constexpr auto BinaryWriter<eT>::is_supported_type() -> bool
{
    return std::is_same_v<eT, uint8_t> || std::is_same_v<eT, float>
           || std::is_same_v<eT, double>;
}

template <typename eT>
constexpr auto BinaryWriter<eT>::dtype_code() -> uint8_t
{
    if constexpr (std::is_same_v<eT, uint8_t>)
    {
        return 1;
    }
    if constexpr (std::is_same_v<eT, float>)
    {
        return 2;
    }
    return 3;
}

template <typename eT>
BinaryWriter<eT>::BinaryWriter(std::string_view file_path)
    : path_(std::string(file_path)), io_buffer_(kDefaultBufferSize)
{
    static_assert(
        is_supported_type(),
        "BinaryWriter only supports uint8_t, float, and double");

    file_ = detail::open_file<std::ofstream>(
        path_, std::ios::binary | std::ios::trunc, io_buffer_);

    write_meta();
}

template <typename eT>
BinaryWriter<eT>::~BinaryWriter() noexcept
{
    try
    {
        finish();
    }
    catch (...)
    {
    }
}

template <typename eT>
auto BinaryWriter<eT>::write(
    const Eigen::Ref<const Eigen::Matrix<eT, Eigen::Dynamic, 1>>& record)
    -> void
{
    if (finished_)
    {
        throw InvalidOperationException(
            std::format("{}: cannot write after finish", path_.string()));
    }

    if (!has_shape_)
    {
        n_rows_ = static_cast<uint64_t>(record.size());
        has_shape_ = true;
    }
    else if (n_rows_ != static_cast<uint64_t>(record.size()))
    {
        throw ArgumentValidationException(
            std::format(
                "{}: inconsistent record size, expected {}, got {}",
                path_.string(),
                n_rows_,
                record.size()));
    }

    const auto bytes_count = static_cast<size_t>(record.size()) * sizeof(eT);
    if (bytes_count > 0)
    {
        file_.write(
            reinterpret_cast<const char*>(record.data()),
            static_cast<std::streamsize>(bytes_count));

        if (!file_.good())
        {
            throw FileWriteException(
                std::format("{}: failed to write record data", path_.string()));
        }
    }

    ++n_cols_;
}

template <typename eT>
auto BinaryWriter<eT>::finish() -> void
{
    if (finished_)
    {
        return;
    }

    file_.seekp(0, std::ios::beg);
    if (!file_.good())
    {
        throw FileWriteException(
            std::format("{}: failed to seek to file header", path_.string()));
    }

    write_meta();

    file_.flush();
    if (!file_.good())
    {
        throw FileWriteException(
            std::format("{}: failed to flush output file", path_.string()));
    }

    finished_ = true;
}

template <typename eT>
auto BinaryWriter<eT>::write_meta() -> void
{
    file_.write(
        reinterpret_cast<const char*>(kMagic.data()),
        static_cast<std::streamsize>(kMagic.size()));
    write_u32_le(kVersion);
    write_u64_le(n_rows_);
    write_u64_le(n_cols_);

    const auto dtype = dtype_code();
    file_.write(reinterpret_cast<const char*>(&dtype), 1);

    if (!file_.good())
    {
        throw FileWriteException(
            std::format("{}: failed to write metadata", path_.string()));
    }
}

template <typename eT>
auto BinaryWriter<eT>::write_u32_le(uint32_t value) -> void
{
    std::array<unsigned char, 4> bytes
        = {static_cast<unsigned char>(value & 0xFFU),
           static_cast<unsigned char>((value >> 8U) & 0xFFU),
           static_cast<unsigned char>((value >> 16U) & 0xFFU),
           static_cast<unsigned char>((value >> 24U) & 0xFFU)};

    file_.write(
        reinterpret_cast<const char*>(bytes.data()),
        static_cast<std::streamsize>(bytes.size()));
}

template <typename eT>
auto BinaryWriter<eT>::write_u64_le(uint64_t value) -> void
{
    std::array<unsigned char, 8> bytes
        = {static_cast<unsigned char>(value & 0xFFU),
           static_cast<unsigned char>((value >> 8U) & 0xFFU),
           static_cast<unsigned char>((value >> 16U) & 0xFFU),
           static_cast<unsigned char>((value >> 24U) & 0xFFU),
           static_cast<unsigned char>((value >> 32U) & 0xFFU),
           static_cast<unsigned char>((value >> 40U) & 0xFFU),
           static_cast<unsigned char>((value >> 48U) & 0xFFU),
           static_cast<unsigned char>((value >> 56U) & 0xFFU)};

    file_.write(
        reinterpret_cast<const char*>(bytes.data()),
        static_cast<std::streamsize>(bytes.size()));
}

}  // namespace gelex::detail

#endif  // GELEX_DATA_BINARY_WRITER_H_
