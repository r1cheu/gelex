// Copyright 2026 RuLei Chen
// SPDX-License-Identifier: Apache-2.0

#ifndef GELEX_DATA_DETAIL_BED_SOURCE_H_
#define GELEX_DATA_DETAIL_BED_SOURCE_H_

#include <cstddef>
#include <cstdint>
#include <filesystem>
#include <fmt/format.h>
#include <span>
#include <string>

#include "gelex/exception.h"
#include "gelex/io/mapped_file.h"

namespace gelex::detail
{

class BedSource
{
   public:
    using value_type = std::uint8_t;
    using size_type = std::size_t;
    using const_pointer_type = const std::uint8_t*;
    using const_span_type = std::span<const std::uint8_t>;

    BedSource(
        const std::filesystem::path& bed_path,
        size_type num_variants,
        size_type num_samples);

    BedSource(const BedSource&) = delete;
    BedSource& operator=(const BedSource&) = delete;
    BedSource(BedSource&&) noexcept = default;
    BedSource& operator=(BedSource&&) noexcept = default;
    ~BedSource() = default;

    [[nodiscard]] auto data() const noexcept -> const_pointer_type
    {
        return payload_;
    }

    [[nodiscard]] auto size() const noexcept -> size_type
    {
        return num_variants_;
    }

    [[nodiscard]] auto empty() const noexcept -> bool
    {
        return num_variants_ == 0;
    }

    [[nodiscard]] auto num_samples() const noexcept -> size_type
    {
        return num_samples_;
    }

    [[nodiscard]] auto stride() const noexcept -> size_type
    {
        return bytes_per_variant_;
    }

    [[nodiscard]] auto size_bytes() const noexcept -> size_type
    {
        return num_variants_ * bytes_per_variant_;
    }

    [[nodiscard]] auto bytes() const noexcept -> const_span_type
    {
        return {payload_, size_bytes()};
    }

    [[nodiscard]] auto operator[](size_type variant_idx) const noexcept
        -> const_span_type
    {
        return {
            payload_ + (variant_idx * bytes_per_variant_),
            bytes_per_variant_,
        };
    }

    [[nodiscard]] auto at(size_type variant_idx) const -> const_span_type
    {
        if (variant_idx >= num_variants_)
        {
            throw GelexException(
                fmt::format(
                    "BED variant index out of range: index={}, size={}",
                    variant_idx,
                    num_variants_));
        }

        return (*this)[variant_idx];
    }

    [[nodiscard]] auto subspan(size_type start, size_type count) const
        -> const_span_type
    {
        if (start > num_variants_ || count > num_variants_ - start)
        {
            throw GelexException(
                fmt::format(
                    "BED variant range out of bounds: start={}, count={}, "
                    "size={}",
                    start,
                    count,
                    num_variants_));
        }

        return {
            payload_ + (start * bytes_per_variant_),
            count * bytes_per_variant_,
        };
    }

   private:
    MappedFile mmap_;
    const std::uint8_t* payload_ = nullptr;
    size_type num_variants_ = 0;
    size_type num_samples_ = 0;
    size_type bytes_per_variant_ = 0;
};

[[nodiscard]] auto open_bed_source(const std::string& bfile_prefix)
    -> BedSource;

}  // namespace gelex::detail

#endif  // GELEX_DATA_DETAIL_BED_SOURCE_H_
