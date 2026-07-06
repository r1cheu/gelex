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

#ifndef GELEX_DATA_BED_H_
#define GELEX_DATA_BED_H_

#include <array>
#include <concepts>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <limits>
#include <ranges>
#include <span>
#include <string>
#include <utility>
#include <vector>

#include <fmt/format.h>
#include <Eigen/Core>

#include "gelex/data/dataframe/index.h"
#include "gelex/data/detail/bed_source.h"
#include "gelex/data/detail/index_projection.h"
#include "gelex/exception.h"

namespace gelex
{

class Bed
{
   public:
    Bed(const Bed&) = delete;
    Bed& operator=(const Bed&) = delete;
    Bed(Bed&&) noexcept = default;
    Bed& operator=(Bed&&) noexcept = default;
    ~Bed() = default;

    template <std::floating_point T>
    [[nodiscard]] auto read() const -> Eigen::MatrixX<T>
    {
        return read<T>(0, num_snps());
    }

    template <std::floating_point T>
    [[nodiscard]] auto read(Eigen::Index first_snp, Eigen::Index last_snp) const
        -> Eigen::MatrixX<T>;

    template <std::floating_point T>
    auto read_into(Eigen::Ref<Eigen::MatrixX<T>> out, Eigen::Index first_snp)
        const -> void;

    template <std::floating_point T>
    [[nodiscard]] auto read_snps(
        std::span<const Eigen::Index> snp_indices) const -> Eigen::MatrixX<T>;

    template <std::floating_point T>
    auto read_snps_into(
        Eigen::Ref<Eigen::MatrixX<T>> out,
        std::span<const Eigen::Index> snp_indices) const -> void;

    template <std::floating_point T>
    [[nodiscard]] auto read_snps(const DataFrameIndex<std::string>& target_snps)
        const -> Eigen::MatrixX<T>;

    template <std::floating_point T>
    auto read_snps_into(
        Eigen::Ref<Eigen::MatrixX<T>> out,
        const DataFrameIndex<std::string>& target_snps) const -> void;

    [[nodiscard]] auto num_samples() const noexcept -> Eigen::Index
    {
        return index_projection_.target_size();
    }

    [[nodiscard]] auto num_snps() const noexcept -> Eigen::Index
    {
        return static_cast<Eigen::Index>(bed_source_.size());
    }

    [[nodiscard]] auto sample_index() const noexcept
        -> const DataFrameIndex<std::string>&
    {
        return sample_index_;
    }

    [[nodiscard]] auto snp_index() const noexcept
        -> const DataFrameIndex<std::string>&
    {
        return snp_index_;
    }

   private:
    Bed(detail::BedSource bed_source,
        detail::IndexProjection index_projection,
        DataFrameIndex<std::string> sample_index,
        DataFrameIndex<std::string> snp_index)
        : bed_source_{std::move(bed_source)},
          index_projection_{std::move(index_projection)},
          sample_index_{std::move(sample_index)},
          snp_index_{std::move(snp_index)}
    {
    }

    detail::BedSource bed_source_;
    detail::IndexProjection index_projection_;
    DataFrameIndex<std::string> sample_index_;
    DataFrameIndex<std::string> snp_index_;

    static constexpr std::size_t BED_LUT_SIZE = 256;

    template <std::floating_point T>
    using BedLutEntry = std::array<T, 4>;

    template <std::floating_point T>
    static consteval auto make_bed_lut_entry(std::uint8_t byte)
        -> BedLutEntry<T>;

    template <std::floating_point T>
    static consteval auto make_bed_decode_lut()
        -> std::array<BedLutEntry<T>, BED_LUT_SIZE>;

    template <std::floating_point T>
    static constexpr auto BED_DECODE_LUT = make_bed_decode_lut<T>();

    template <std::floating_point T>
    static auto decode_variant(const std::uint8_t* data_ptr, std::span<T> out)
        -> void;

    friend auto open_bed(const std::string& bfile_prefix) -> Bed;

    friend auto open_bed(
        const std::string& bfile_prefix,
        const DataFrameIndex<std::string>& target_index) -> Bed;
};

[[nodiscard]] auto open_bed(const std::string& bfile_prefix) -> Bed;

[[nodiscard]] auto open_bed(
    const std::string& bfile_prefix,
    const DataFrameIndex<std::string>& target_index) -> Bed;

template <std::floating_point T>
consteval auto Bed::make_bed_lut_entry(std::uint8_t byte) -> BedLutEntry<T>
{
    constexpr BedLutEntry<T> dosage_map{
        T{2},
        std::numeric_limits<T>::quiet_NaN(),
        T{1},
        T{0},
    };

    BedLutEntry<T> entry{};

    for (int i = 0; i < 4; ++i)
    {
        const auto code = static_cast<std::size_t>((byte >> (2 * i)) & 0x03);
        entry[static_cast<std::size_t>(i)] = dosage_map[code];
    }

    return entry;
}

template <std::floating_point T>
consteval auto Bed::make_bed_decode_lut()
    -> std::array<BedLutEntry<T>, BED_LUT_SIZE>
{
    std::array<BedLutEntry<T>, BED_LUT_SIZE> table{};

    for (std::size_t i = 0; i < BED_LUT_SIZE; ++i)
    {
        table[i] = make_bed_lut_entry<T>(static_cast<std::uint8_t>(i));
    }

    return table;
}

template <std::floating_point T>
auto Bed::decode_variant(const std::uint8_t* data_ptr, std::span<T> out) -> void
{
    const auto num_samples = static_cast<Eigen::Index>(out.size());

    const Eigen::Index full_bytes = num_samples / 4;
    const Eigen::Index tail = num_samples % 4;

    for (Eigen::Index i = 0; i < full_bytes; ++i)
    {
        const auto& vals = BED_DECODE_LUT<T>[data_ptr[i]];

        std::memcpy(out.data() + (i * 4), vals.data(), 4 * sizeof(T));
    }

    if (tail != 0)
    {
        const auto& vals = BED_DECODE_LUT<T>[data_ptr[full_bytes]];
        const Eigen::Index base = full_bytes * 4;

        for (Eigen::Index k = 0; k < tail; ++k)
        {
            out[static_cast<std::size_t>(base + k)]
                = vals[static_cast<std::size_t>(k)];
        }
    }
}

template <std::floating_point T>
auto Bed::read(Eigen::Index first_snp, Eigen::Index last_snp) const
    -> Eigen::MatrixX<T>
{
    if (first_snp < 0 || last_snp < first_snp || last_snp > num_snps())
    {
        throw GelexException(
            fmt::format(
                "invalid BED SNP range: [{}, {}). Total SNPs: {}",
                first_snp,
                last_snp,
                num_snps()));
    }

    Eigen::MatrixX<T> out(num_samples(), last_snp - first_snp);
    read_into<T>(out, first_snp);
    return out;
}

template <std::floating_point T>
auto Bed::read_into(Eigen::Ref<Eigen::MatrixX<T>> out, Eigen::Index first_snp)
    const -> void
{
    if (out.rows() != num_samples())
    {
        throw GelexException(
            fmt::format(
                "BED target buffer row mismatch: expected {}, got {}",
                num_samples(),
                out.rows()));
    }

    const auto last_snp = first_snp + out.cols();
    if (first_snp < 0 || last_snp < first_snp || last_snp > num_snps())
    {
        throw GelexException(
            fmt::format(
                "invalid BED SNP range: [{}, {}). Total SNPs: {}",
                first_snp,
                last_snp,
                num_snps()));
    }

    const auto target_to_source = index_projection_.target_to_source();

#pragma omp parallel
    {
        Eigen::VectorX<T> source_values(index_projection_.source_size());

#pragma omp for schedule(static)
        for (Eigen::Index target_col = 0; target_col < out.cols(); ++target_col)
        {
            const auto source_variant = first_snp + target_col;
            const auto source_bytes
                = bed_source_[static_cast<std::size_t>(source_variant)];

            decode_variant(
                source_bytes.data(),
                std::span<T>{
                    source_values.data(),
                    static_cast<std::size_t>(source_values.size())});

            for (const auto [target_row, source_row] :
                 std::views::enumerate(target_to_source))
            {
                out(static_cast<Eigen::Index>(target_row), target_col)
                    = source_values[source_row];
            }
        }
    }
}

template <std::floating_point T>
auto Bed::read_snps(std::span<const Eigen::Index> snp_indices) const
    -> Eigen::MatrixX<T>
{
    Eigen::MatrixX<T> out(
        num_samples(), static_cast<Eigen::Index>(snp_indices.size()));
    read_snps_into<T>(out, snp_indices);
    return out;
}

template <std::floating_point T>
auto Bed::read_snps_into(
    Eigen::Ref<Eigen::MatrixX<T>> out,
    std::span<const Eigen::Index> snp_indices) const -> void
{
    if (out.rows() != num_samples())
    {
        throw GelexException(
            fmt::format(
                "BED target buffer row mismatch: expected {}, got {}",
                num_samples(),
                out.rows()));
    }

    if (out.cols() != static_cast<Eigen::Index>(snp_indices.size()))
    {
        throw GelexException(
            fmt::format(
                "BED target buffer column mismatch: expected {}, got {}",
                snp_indices.size(),
                out.cols()));
    }

    for (const auto source_variant : snp_indices)
    {
        if (source_variant < 0 || source_variant >= num_snps())
        {
            throw GelexException(
                fmt::format(
                    "BED SNP index out of range: index={}, size={}",
                    source_variant,
                    num_snps()));
        }
    }

    const auto target_to_source = index_projection_.target_to_source();

#pragma omp parallel
    {
        Eigen::VectorX<T> source_values(index_projection_.source_size());

#pragma omp for schedule(static)
        for (Eigen::Index target_col = 0; target_col < out.cols(); ++target_col)
        {
            const auto source_variant
                = snp_indices[static_cast<std::size_t>(target_col)];
            const auto source_bytes
                = bed_source_[static_cast<std::size_t>(source_variant)];

            decode_variant(
                source_bytes.data(),
                std::span<T>{
                    source_values.data(),
                    static_cast<std::size_t>(source_values.size())});

            for (const auto [target_row, source_row] :
                 std::views::enumerate(target_to_source))
            {
                out(static_cast<Eigen::Index>(target_row), target_col)
                    = source_values[source_row];
            }
        }
    }
}

template <std::floating_point T>
auto Bed::read_snps(const DataFrameIndex<std::string>& target_snps) const
    -> Eigen::MatrixX<T>
{
    Eigen::MatrixX<T> out(
        num_samples(), static_cast<Eigen::Index>(target_snps.size()));
    read_snps_into<T>(out, target_snps);
    return out;
}

template <std::floating_point T>
auto Bed::read_snps_into(
    Eigen::Ref<Eigen::MatrixX<T>> out,
    const DataFrameIndex<std::string>& target_snps) const -> void
{
    if (out.cols() != static_cast<Eigen::Index>(target_snps.size()))
    {
        throw GelexException(
            fmt::format(
                "BED target buffer column mismatch: expected {}, got {}",
                target_snps.size(),
                out.cols()));
    }

    std::vector<Eigen::Index> snp_indices;
    snp_indices.reserve(target_snps.size());

    for (const auto& snp : target_snps.keys())
    {
        snp_indices.push_back(static_cast<Eigen::Index>(snp_index_.at(snp)));
    }

    read_snps_into<T>(out, snp_indices);
}

}  // namespace gelex

#endif  // GELEX_DATA_BED_H_
