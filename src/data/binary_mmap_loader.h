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

#ifndef GELEX_DATA_BINARY_MMAP_LOADER_H_
#define GELEX_DATA_BINARY_MMAP_LOADER_H_

#include <array>
#include <cstddef>
#include <cstdint>
#include <filesystem>
#include <format>
#include <fstream>
#include <limits>
#include <string>
#include <string_view>
#include <system_error>
#include <type_traits>

#include <mio.h>

#include <Eigen/Core>

#include "gelex/exception.h"

namespace gelex::detail
{

template <typename eT>
class BinaryMmapLoader
{
   public:
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
    static constexpr size_t kMetaSize = 32;

    using MatrixType = Eigen::Matrix<eT, Eigen::Dynamic, Eigen::Dynamic>;
    using MapType = Eigen::Map<const MatrixType, Eigen::Unaligned>;

    explicit BinaryMmapLoader(std::string_view file_path);

    BinaryMmapLoader(const BinaryMmapLoader&) = delete;
    BinaryMmapLoader(BinaryMmapLoader&&) noexcept = default;
    auto operator=(const BinaryMmapLoader&) -> BinaryMmapLoader& = delete;
    auto operator=(BinaryMmapLoader&&) noexcept -> BinaryMmapLoader& = default;
    ~BinaryMmapLoader() = default;

    [[nodiscard]] auto matrix() const noexcept -> const MapType&;
    [[nodiscard]] auto load_copy() const -> MatrixType;

   private:
    static constexpr auto is_supported_type() -> bool;
    static constexpr auto dtype_code() -> uint8_t;

    static auto decode_u32_le(const std::byte* data) -> uint32_t;
    static auto decode_u64_le(const std::byte* data) -> uint64_t;

    auto checked_mul(size_t a, size_t b, std::string_view ctx) const -> size_t;
    auto checked_add(size_t a, size_t b, std::string_view ctx) const -> size_t;

    auto parse_and_validate_meta() -> void;

    std::filesystem::path path_;
    mio::mmap_source mmap_;
    MapType matrix_view_;
    uint64_t n_rows_ = 0;
    uint64_t n_cols_ = 0;
    size_t payload_bytes_ = 0;
};

template <typename eT>
constexpr auto BinaryMmapLoader<eT>::is_supported_type() -> bool
{
    return std::is_same_v<eT, uint8_t> || std::is_same_v<eT, float>
           || std::is_same_v<eT, double>;
}

template <typename eT>
constexpr auto BinaryMmapLoader<eT>::dtype_code() -> uint8_t
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
BinaryMmapLoader<eT>::BinaryMmapLoader(std::string_view file_path)
    : path_(std::string(file_path)), matrix_view_(nullptr, 0, 0)
{
    static_assert(
        is_supported_type(),
        "BinaryMmapLoader only supports uint8_t, float, and double");

    parse_and_validate_meta();

    if (payload_bytes_ > 0)
    {
        std::error_code ec;
        mmap_.map(path_.string(), kMetaSize, payload_bytes_, ec);
        if (ec)
        {
            throw FileOpenException(
                std::format(
                    "{}: failed to mmap payload: {}",
                    path_.string(),
                    ec.message()));
        }
    }

    const auto* payload_ptr = reinterpret_cast<const eT*>(mmap_.data());
    new (&matrix_view_) MapType(
        payload_ptr,
        static_cast<Eigen::Index>(n_rows_),
        static_cast<Eigen::Index>(n_cols_));
}

template <typename eT>
auto BinaryMmapLoader<eT>::matrix() const noexcept -> const MapType&
{
    return matrix_view_;
}

template <typename eT>
auto BinaryMmapLoader<eT>::load_copy() const -> MatrixType
{
    return MatrixType(matrix_view_);
}

template <typename eT>
auto BinaryMmapLoader<eT>::decode_u32_le(const std::byte* data) -> uint32_t
{
    return static_cast<uint32_t>(static_cast<uint8_t>(data[0]))
           | (static_cast<uint32_t>(static_cast<uint8_t>(data[1])) << 8U)
           | (static_cast<uint32_t>(static_cast<uint8_t>(data[2])) << 16U)
           | (static_cast<uint32_t>(static_cast<uint8_t>(data[3])) << 24U);
}

template <typename eT>
auto BinaryMmapLoader<eT>::decode_u64_le(const std::byte* data) -> uint64_t
{
    uint64_t value = 0;
    for (size_t i = 0; i < 8; ++i)
    {
        value |= static_cast<uint64_t>(static_cast<uint8_t>(data[i]))
                 << (i * 8U);
    }
    return value;
}

template <typename eT>
auto BinaryMmapLoader<eT>::checked_mul(size_t a, size_t b, std::string_view ctx)
    const -> size_t
{
    if (a == 0 || b == 0)
    {
        return 0;
    }

    if (a > (std::numeric_limits<size_t>::max() / b))
    {
        throw FileFormatException(
            std::format(
                "{}: size overflow while computing {}", path_.string(), ctx));
    }

    return a * b;
}

template <typename eT>
auto BinaryMmapLoader<eT>::checked_add(size_t a, size_t b, std::string_view ctx)
    const -> size_t
{
    if (a > (std::numeric_limits<size_t>::max() - b))
    {
        throw FileFormatException(
            std::format(
                "{}: size overflow while computing {}", path_.string(), ctx));
    }

    return a + b;
}

template <typename eT>
auto BinaryMmapLoader<eT>::parse_and_validate_meta() -> void
{
    std::ifstream file(path_, std::ios::binary);
    if (!file.is_open())
    {
        if (!std::filesystem::exists(path_))
        {
            throw FileNotFoundException(
                std::format("{}: not found", path_.string()));
        }

        throw FileOpenException(
            std::format("{}: failed to open file", path_.string()));
    }

    std::array<std::byte, kMetaSize> header{};
    file.read(
        reinterpret_cast<char*>(header.data()),
        static_cast<std::streamsize>(kMetaSize));
    if (!file || file.gcount() != static_cast<std::streamsize>(kMetaSize))
    {
        throw FileFormatException(
            std::format(
                "{}: file too small for binary header", path_.string()));
    }

    const auto* raw = header.data();
    for (size_t i = 0; i < kMagic.size(); ++i)
    {
        if (raw[i] != kMagic[i])
        {
            throw FileFormatException(
                std::format("{}: invalid file magic", path_.string()));
        }
    }

    const auto version = decode_u32_le(raw + 8);
    if (version != kVersion)
    {
        throw FileFormatException(
            std::format(
                "{}: unsupported version {}, expected {}",
                path_.string(),
                version,
                kVersion));
    }

    n_rows_ = decode_u64_le(raw + 12);
    n_cols_ = decode_u64_le(raw + 20);
    const auto stored_dtype = static_cast<uint8_t>(raw[28]);

    if (stored_dtype < 1 || stored_dtype > 3)
    {
        throw FileFormatException(
            std::format(
                "{}: invalid dtype code {}", path_.string(), stored_dtype));
    }

    if (stored_dtype != dtype_code())
    {
        throw ArgumentValidationException(
            std::format(
                "{}: dtype mismatch, file={}, requested={}",
                path_.string(),
                stored_dtype,
                dtype_code()));
    }

    const auto rows = static_cast<size_t>(n_rows_);
    const auto cols = static_cast<size_t>(n_cols_);
    const auto total_elements = checked_mul(rows, cols, "number of elements");
    payload_bytes_ = checked_mul(total_elements, sizeof(eT), "payload bytes");
    const auto expected_size
        = checked_add(kMetaSize, payload_bytes_, "file bytes");

    std::error_code file_size_ec;
    const auto actual_size = std::filesystem::file_size(path_, file_size_ec);
    if (file_size_ec)
    {
        throw FileOpenException(
            std::format("{}: failed to query file size", path_.string()));
    }

    if (actual_size != expected_size)
    {
        throw FileFormatException(
            std::format(
                "{}: file size mismatch, expected {} bytes, got {} bytes",
                path_.string(),
                expected_size,
                actual_size));
    }
}

}  // namespace gelex::detail

#endif  // GELEX_DATA_BINARY_MMAP_LOADER_H_
