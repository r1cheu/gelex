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

#include <Eigen/Core>
#include <array>
#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>
#include <cblas.h>
#include <chrono>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <nanobench.h>
#include <random>
#include <string>
#include <string_view>
#include <vector>

#include "bayes/detail/genotype_kernel.h"

#if defined(__x86_64__) || defined(_M_X64)
#include <immintrin.h>
#define GELEX_BENCHMARK_X86_64 1
#else
#define GELEX_BENCHMARK_X86_64 0
#endif

#if GELEX_BENCHMARK_X86_64 && (defined(__GNUC__) || defined(__clang__))
#define GELEX_AVX2_TARGET __attribute__((target("avx2,fma")))
#define GELEX_NOINLINE __attribute__((noinline))
#else
#define GELEX_AVX2_TARGET
#define GELEX_NOINLINE
#endif

namespace
{

using LocusLut = std::array<double, 4>;

consteval auto make_bed_unpack_table() -> std::array<std::uint32_t, 256>
{
    std::array<std::uint32_t, 256> table{};
    for (std::uint32_t packed = 0; packed < table.size(); ++packed)
    {
        table[packed] = (packed & 0x03U) | (((packed >> 2U) & 0x03U) << 8U)
                        | (((packed >> 4U) & 0x03U) << 16U)
                        | (((packed >> 6U) & 0x03U) << 24U);
    }
    return table;
}

inline constexpr auto BED_UNPACK_TABLE = make_bed_unpack_table();

struct Shape
{
    Eigen::Index samples;
    Eigen::Index markers;
};

struct CompactGenotypeFixture
{
    CompactGenotypeFixture(Eigen::Index samples, Eigen::Index markers)
        : samples(samples),
          markers(markers),
          bed_stride(static_cast<std::size_t>((samples + 3) / 4)),
          dense(samples, markers),
          uint8_codes(
              static_cast<std::size_t>(samples)
              * static_cast<std::size_t>(markers)),
          bed_codes(bed_stride * static_cast<std::size_t>(markers), 0),
          luts(static_cast<std::size_t>(markers)),
          alpha(static_cast<std::size_t>(markers)),
          initial_residual(samples)
    {
        std::mt19937_64 rng{42};
        std::uniform_real_distribution<double> uniform{0.0, 1.0};
        std::normal_distribution<double> normal{0.0, 1.0};

        for (Eigen::Index marker = 0; marker < markers; ++marker)
        {
            const double phase = static_cast<double>(
                (static_cast<std::uint64_t>(marker) * 104729ULL) % 1000ULL);
            const double allele_frequency = 0.05 + (0.45 * phase / 999.0);
            const double mean = 2.0 * allele_frequency;
            const double inverse_sd
                = 1.0
                  / std::sqrt(
                      2.0 * allele_frequency * (1.0 - allele_frequency));
            auto& lut = luts[static_cast<std::size_t>(marker)];
            lut
                = {(2.0 - mean) * inverse_sd,
                   0.0,
                   (1.0 - mean) * inverse_sd,
                   -mean * inverse_sd};
            alpha[static_cast<std::size_t>(marker)]
                = 1e-4 * static_cast<double>((marker % 7) - 3);

            auto* unpacked = uint8_column(marker);
            auto* packed = bed_column(marker);
            const double p_hom_a1 = allele_frequency * allele_frequency;
            const double p_het
                = 2.0 * allele_frequency * (1.0 - allele_frequency);
            for (Eigen::Index sample = 0; sample < samples; ++sample)
            {
                const double draw = uniform(rng);
                std::uint8_t code = 3;
                if (draw < 0.01)
                {
                    code = 1;
                }
                else
                {
                    const double genotype_draw = (draw - 0.01) / 0.99;
                    if (genotype_draw < p_hom_a1)
                    {
                        code = 0;
                    }
                    else if (genotype_draw < p_hom_a1 + p_het)
                    {
                        code = 2;
                    }
                }

                unpacked[sample] = code;
                const auto byte_index = static_cast<std::size_t>(sample / 4);
                const auto shift = static_cast<unsigned>(2 * (sample % 4));
                packed[byte_index] |= static_cast<std::uint8_t>(code << shift);
                dense(sample, marker) = lut[code];
            }
        }

        for (Eigen::Index sample = 0; sample < samples; ++sample)
        {
            initial_residual(sample) = normal(rng);
        }
    }

    [[nodiscard]] auto uint8_column(Eigen::Index marker) noexcept
        -> std::uint8_t*
    {
        return uint8_codes.data()
               + (static_cast<std::size_t>(marker)
                  * static_cast<std::size_t>(samples));
    }

    [[nodiscard]] auto uint8_column(Eigen::Index marker) const noexcept
        -> const std::uint8_t*
    {
        return uint8_codes.data()
               + (static_cast<std::size_t>(marker)
                  * static_cast<std::size_t>(samples));
    }

    [[nodiscard]] auto bed_column(Eigen::Index marker) noexcept -> std::uint8_t*
    {
        return bed_codes.data()
               + (static_cast<std::size_t>(marker) * bed_stride);
    }

    [[nodiscard]] auto bed_column(Eigen::Index marker) const noexcept
        -> const std::uint8_t*
    {
        return bed_codes.data()
               + (static_cast<std::size_t>(marker) * bed_stride);
    }

    Eigen::Index samples;
    Eigen::Index markers;
    std::size_t bed_stride;
    Eigen::MatrixXd dense;
    std::vector<std::uint8_t> uint8_codes;
    std::vector<std::uint8_t> bed_codes;
    std::vector<LocusLut> luts;
    std::vector<double> alpha;
    Eigen::VectorXd initial_residual;
};

GELEX_NOINLINE auto dot_bed_scalar(
    const std::uint8_t* packed,
    const double* lut,
    const double* residual,
    Eigen::Index size) noexcept -> double
{
    double sum = 0.0;
    for (Eigen::Index i = 0; i < size; ++i)
    {
        const auto byte = packed[static_cast<std::size_t>(i / 4)];
        const auto shift = static_cast<unsigned>(2 * (i % 4));
        sum += lut[(byte >> shift) & 0x03U] * residual[i];
    }
    return sum;
}

GELEX_NOINLINE auto axpy_bed_scalar(
    const std::uint8_t* packed,
    const double* lut,
    double scale,
    double* residual,
    Eigen::Index size) noexcept -> void
{
    for (Eigen::Index i = 0; i < size; ++i)
    {
        const auto byte = packed[static_cast<std::size_t>(i / 4)];
        const auto shift = static_cast<unsigned>(2 * (i % 4));
        residual[i] += scale * lut[(byte >> shift) & 0x03U];
    }
}

GELEX_NOINLINE auto materialize_dot_uint8_scalar(
    const std::uint8_t* codes,
    const double* lut,
    const double* residual,
    double* scratch,
    Eigen::Index size) noexcept -> double
{
    double sum = 0.0;
    for (Eigen::Index i = 0; i < size; ++i)
    {
        const double value = lut[codes[i]];
        scratch[i] = value;
        sum += value * residual[i];
    }
    return sum;
}

GELEX_NOINLINE auto materialize_dot_bed_scalar(
    const std::uint8_t* packed,
    const double* lut,
    const double* residual,
    double* scratch,
    Eigen::Index size) noexcept -> double
{
    double sum = 0.0;
    for (Eigen::Index i = 0; i < size; ++i)
    {
        const auto byte = packed[static_cast<std::size_t>(i / 4)];
        const auto shift = static_cast<unsigned>(2 * (i % 4));
        const double value = lut[(byte >> shift) & 0x03U];
        scratch[i] = value;
        sum += value * residual[i];
    }
    return sum;
}

#if GELEX_BENCHMARK_X86_64 && (defined(__GNUC__) || defined(__clang__))

GELEX_AVX2_TARGET inline auto lookup4(
    std::uint32_t four_codes,
    const double* lut) noexcept -> __m256d
{
    const __m128i code_bytes
        = _mm_cvtsi32_si128(static_cast<std::int32_t>(four_codes));
    const __m128i codes = _mm_cvtepu8_epi32(code_bytes);
    const __m128i even_indices = _mm_slli_epi32(codes, 1);
    const __m128i odd_indices = _mm_add_epi32(even_indices, _mm_set1_epi32(1));
    const __m256i indices = _mm256_set_m128i(
        _mm_unpackhi_epi32(even_indices, odd_indices),
        _mm_unpacklo_epi32(even_indices, odd_indices));
    const __m256i lut_parts
        = _mm256_loadu_si256(reinterpret_cast<const __m256i*>(lut));
    return _mm256_castsi256_pd(_mm256_permutevar8x32_epi32(lut_parts, indices));
}

GELEX_AVX2_TARGET inline auto load_uint8_values(
    const std::uint8_t* codes,
    const double* lut) noexcept -> __m256d
{
    std::uint32_t four_codes = 0;
    std::memcpy(&four_codes, codes, sizeof(four_codes));
    return lookup4(four_codes, lut);
}

GELEX_AVX2_TARGET inline auto load_bed_values(
    std::uint8_t packed,
    const double* lut) noexcept -> __m256d
{
    return lookup4(BED_UNPACK_TABLE[packed], lut);
}

GELEX_AVX2_TARGET inline auto unpack_bed8(const std::uint8_t* packed) noexcept
    -> std::uint64_t
{
    return static_cast<std::uint64_t>(BED_UNPACK_TABLE[packed[0]])
           | (static_cast<std::uint64_t>(BED_UNPACK_TABLE[packed[1]]) << 32U);
}

GELEX_AVX2_TARGET inline auto horizontal_sum(__m256d values) noexcept -> double
{
    const __m128d low = _mm256_castpd256_pd128(values);
    const __m128d high = _mm256_extractf128_pd(values, 1);
    const __m128d pair = _mm_add_pd(low, high);
    const __m128d swapped = _mm_shuffle_pd(pair, pair, 0x01);
    return _mm_cvtsd_f64(_mm_add_sd(pair, swapped));
}

GELEX_AVX2_TARGET GELEX_NOINLINE auto dot_bed_avx2(
    const std::uint8_t* packed,
    const double* lut,
    const double* residual,
    Eigen::Index size) noexcept -> double
{
    __m256d sum0 = _mm256_setzero_pd();
    __m256d sum1 = _mm256_setzero_pd();
    __m256d sum2 = _mm256_setzero_pd();
    __m256d sum3 = _mm256_setzero_pd();
    Eigen::Index i = 0;
    std::size_t byte = 0;
    for (; i + 16 <= size; i += 16, byte += 4)
    {
        const std::uint64_t codes0 = unpack_bed8(packed + byte);
        const std::uint64_t codes1 = unpack_bed8(packed + byte + 2);
        sum0 = _mm256_fmadd_pd(
            lookup4(static_cast<std::uint32_t>(codes0), lut),
            _mm256_loadu_pd(residual + i),
            sum0);
        sum1 = _mm256_fmadd_pd(
            lookup4(static_cast<std::uint32_t>(codes0 >> 32U), lut),
            _mm256_loadu_pd(residual + i + 4),
            sum1);
        sum2 = _mm256_fmadd_pd(
            lookup4(static_cast<std::uint32_t>(codes1), lut),
            _mm256_loadu_pd(residual + i + 8),
            sum2);
        sum3 = _mm256_fmadd_pd(
            lookup4(static_cast<std::uint32_t>(codes1 >> 32U), lut),
            _mm256_loadu_pd(residual + i + 12),
            sum3);
    }
    __m256d sum
        = _mm256_add_pd(_mm256_add_pd(sum0, sum1), _mm256_add_pd(sum2, sum3));
    for (; i + 4 <= size; i += 4, ++byte)
    {
        sum = _mm256_fmadd_pd(
            load_bed_values(packed[byte], lut),
            _mm256_loadu_pd(residual + i),
            sum);
    }
    double result = horizontal_sum(sum);
    for (; i < size; ++i)
    {
        const auto packed_byte = packed[static_cast<std::size_t>(i / 4)];
        const auto shift = static_cast<unsigned>(2 * (i % 4));
        result += lut[(packed_byte >> shift) & 0x03U] * residual[i];
    }
    return result;
}

GELEX_AVX2_TARGET GELEX_NOINLINE auto axpy_bed_avx2(
    const std::uint8_t* packed,
    const double* lut,
    double scale,
    double* residual,
    Eigen::Index size) noexcept -> void
{
    const __m256d scale4 = _mm256_set1_pd(scale);
    Eigen::Index i = 0;
    std::size_t byte = 0;
    for (; i + 8 <= size; i += 8, byte += 2)
    {
        const std::uint64_t codes = unpack_bed8(packed + byte);
        const __m256d updated0 = _mm256_fmadd_pd(
            scale4,
            lookup4(static_cast<std::uint32_t>(codes), lut),
            _mm256_loadu_pd(residual + i));
        const __m256d updated1 = _mm256_fmadd_pd(
            scale4,
            lookup4(static_cast<std::uint32_t>(codes >> 32U), lut),
            _mm256_loadu_pd(residual + i + 4));
        _mm256_storeu_pd(residual + i, updated0);
        _mm256_storeu_pd(residual + i + 4, updated1);
    }
    for (; i + 4 <= size; i += 4, ++byte)
    {
        const __m256d updated = _mm256_fmadd_pd(
            scale4,
            load_bed_values(packed[byte], lut),
            _mm256_loadu_pd(residual + i));
        _mm256_storeu_pd(residual + i, updated);
    }
    for (; i < size; ++i)
    {
        const auto packed_byte = packed[static_cast<std::size_t>(i / 4)];
        const auto shift = static_cast<unsigned>(2 * (i % 4));
        residual[i] += scale * lut[(packed_byte >> shift) & 0x03U];
    }
}

GELEX_AVX2_TARGET GELEX_NOINLINE auto materialize_dot_uint8_avx2(
    const std::uint8_t* codes,
    const double* lut,
    const double* residual,
    double* scratch,
    Eigen::Index size) noexcept -> double
{
    __m256d sum0 = _mm256_setzero_pd();
    __m256d sum1 = _mm256_setzero_pd();
    __m256d sum2 = _mm256_setzero_pd();
    __m256d sum3 = _mm256_setzero_pd();
    Eigen::Index i = 0;
    for (; i + 16 <= size; i += 16)
    {
        const __m256d value0 = load_uint8_values(codes + i, lut);
        const __m256d value1 = load_uint8_values(codes + i + 4, lut);
        const __m256d value2 = load_uint8_values(codes + i + 8, lut);
        const __m256d value3 = load_uint8_values(codes + i + 12, lut);
        _mm256_storeu_pd(scratch + i, value0);
        _mm256_storeu_pd(scratch + i + 4, value1);
        _mm256_storeu_pd(scratch + i + 8, value2);
        _mm256_storeu_pd(scratch + i + 12, value3);
        sum0 = _mm256_fmadd_pd(value0, _mm256_loadu_pd(residual + i), sum0);
        sum1 = _mm256_fmadd_pd(value1, _mm256_loadu_pd(residual + i + 4), sum1);
        sum2 = _mm256_fmadd_pd(value2, _mm256_loadu_pd(residual + i + 8), sum2);
        sum3
            = _mm256_fmadd_pd(value3, _mm256_loadu_pd(residual + i + 12), sum3);
    }
    __m256d sum
        = _mm256_add_pd(_mm256_add_pd(sum0, sum1), _mm256_add_pd(sum2, sum3));
    for (; i + 4 <= size; i += 4)
    {
        const __m256d values = load_uint8_values(codes + i, lut);
        _mm256_storeu_pd(scratch + i, values);
        sum = _mm256_fmadd_pd(values, _mm256_loadu_pd(residual + i), sum);
    }
    double result = horizontal_sum(sum);
    for (; i < size; ++i)
    {
        const double value = lut[codes[i]];
        scratch[i] = value;
        result += value * residual[i];
    }
    return result;
}

GELEX_AVX2_TARGET GELEX_NOINLINE auto materialize_dot_bed_avx2(
    const std::uint8_t* packed,
    const double* lut,
    const double* residual,
    double* scratch,
    Eigen::Index size) noexcept -> double
{
    __m256d sum0 = _mm256_setzero_pd();
    __m256d sum1 = _mm256_setzero_pd();
    __m256d sum2 = _mm256_setzero_pd();
    __m256d sum3 = _mm256_setzero_pd();
    Eigen::Index i = 0;
    std::size_t byte = 0;
    for (; i + 16 <= size; i += 16, byte += 4)
    {
        const std::uint64_t codes0 = unpack_bed8(packed + byte);
        const std::uint64_t codes1 = unpack_bed8(packed + byte + 2);
        const __m256d value0 = lookup4(static_cast<std::uint32_t>(codes0), lut);
        const __m256d value1
            = lookup4(static_cast<std::uint32_t>(codes0 >> 32U), lut);
        const __m256d value2 = lookup4(static_cast<std::uint32_t>(codes1), lut);
        const __m256d value3
            = lookup4(static_cast<std::uint32_t>(codes1 >> 32U), lut);
        _mm256_storeu_pd(scratch + i, value0);
        _mm256_storeu_pd(scratch + i + 4, value1);
        _mm256_storeu_pd(scratch + i + 8, value2);
        _mm256_storeu_pd(scratch + i + 12, value3);
        sum0 = _mm256_fmadd_pd(value0, _mm256_loadu_pd(residual + i), sum0);
        sum1 = _mm256_fmadd_pd(value1, _mm256_loadu_pd(residual + i + 4), sum1);
        sum2 = _mm256_fmadd_pd(value2, _mm256_loadu_pd(residual + i + 8), sum2);
        sum3
            = _mm256_fmadd_pd(value3, _mm256_loadu_pd(residual + i + 12), sum3);
    }
    __m256d sum
        = _mm256_add_pd(_mm256_add_pd(sum0, sum1), _mm256_add_pd(sum2, sum3));
    for (; i + 4 <= size; i += 4, ++byte)
    {
        const __m256d values = load_bed_values(packed[byte], lut);
        _mm256_storeu_pd(scratch + i, values);
        sum = _mm256_fmadd_pd(values, _mm256_loadu_pd(residual + i), sum);
    }
    double result = horizontal_sum(sum);
    for (; i < size; ++i)
    {
        const auto packed_byte = packed[static_cast<std::size_t>(i / 4)];
        const auto shift = static_cast<unsigned>(2 * (i % 4));
        const double value = lut[(packed_byte >> shift) & 0x03U];
        scratch[i] = value;
        result += value * residual[i];
    }
    return result;
}

#endif

struct DenseOps
{
    const CompactGenotypeFixture& fixture;

    [[nodiscard]] auto dot(Eigen::Index marker, const double* residual) const
        -> double
    {
        return fixture.dense.col(marker).dot(
            Eigen::Map<const Eigen::VectorXd>(residual, fixture.samples));
    }

    auto axpy(Eigen::Index marker, double scale, double* residual) const -> void
    {
        Eigen::Map<Eigen::VectorXd>(residual, fixture.samples).noalias()
            += scale * fixture.dense.col(marker);
    }
};

struct BlasDenseOps
{
    const CompactGenotypeFixture& fixture;

    [[nodiscard]] auto dot(Eigen::Index marker, const double* residual) const
        -> double
    {
        return cblas_ddot(
            static_cast<blasint>(fixture.samples),
            fixture.dense.col(marker).data(),
            1,
            residual,
            1);
    }

    auto axpy(Eigen::Index marker, double scale, double* residual) const -> void
    {
        cblas_daxpy(
            static_cast<blasint>(fixture.samples),
            scale,
            fixture.dense.col(marker).data(),
            1,
            residual,
            1);
    }
};

struct Uint8ScalarOps
{
    const CompactGenotypeFixture& fixture;

    [[nodiscard]] auto dot(Eigen::Index marker, const double* residual) const
        -> double
    {
        return gelex::bayes::detail::dot_scalar(
            fixture.uint8_column(marker),
            fixture.luts[static_cast<std::size_t>(marker)].data(),
            residual,
            static_cast<std::size_t>(fixture.samples));
    }

    auto axpy(Eigen::Index marker, double scale, double* residual) const -> void
    {
        gelex::bayes::detail::axpy_scalar(
            fixture.uint8_column(marker),
            fixture.luts[static_cast<std::size_t>(marker)].data(),
            scale,
            residual,
            static_cast<std::size_t>(fixture.samples));
    }
};

struct BedScalarOps
{
    const CompactGenotypeFixture& fixture;

    [[nodiscard]] auto dot(Eigen::Index marker, const double* residual) const
        -> double
    {
        return dot_bed_scalar(
            fixture.bed_column(marker),
            fixture.luts[static_cast<std::size_t>(marker)].data(),
            residual,
            fixture.samples);
    }

    auto axpy(Eigen::Index marker, double scale, double* residual) const -> void
    {
        axpy_bed_scalar(
            fixture.bed_column(marker),
            fixture.luts[static_cast<std::size_t>(marker)].data(),
            scale,
            residual,
            fixture.samples);
    }
};

#if GELEX_BENCHMARK_X86_64 && (defined(__GNUC__) || defined(__clang__))

struct Uint8Avx2Ops
{
    const CompactGenotypeFixture& fixture;

    [[nodiscard]] auto dot(Eigen::Index marker, const double* residual) const
        -> double
    {
        return gelex::bayes::detail::dot_avx2(
            fixture.uint8_column(marker),
            fixture.luts[static_cast<std::size_t>(marker)].data(),
            residual,
            static_cast<std::size_t>(fixture.samples));
    }

    auto axpy(Eigen::Index marker, double scale, double* residual) const -> void
    {
        gelex::bayes::detail::axpy_avx2(
            fixture.uint8_column(marker),
            fixture.luts[static_cast<std::size_t>(marker)].data(),
            scale,
            residual,
            static_cast<std::size_t>(fixture.samples));
    }
};

struct Uint8Avx512Ops
{
    const CompactGenotypeFixture& fixture;

    [[nodiscard]] auto dot(Eigen::Index marker, const double* residual) const
        -> double
    {
        return gelex::bayes::detail::dot_avx512(
            fixture.uint8_column(marker),
            fixture.luts[static_cast<std::size_t>(marker)].data(),
            residual,
            static_cast<std::size_t>(fixture.samples));
    }

    auto axpy(Eigen::Index marker, double scale, double* residual) const -> void
    {
        gelex::bayes::detail::axpy_avx512(
            fixture.uint8_column(marker),
            fixture.luts[static_cast<std::size_t>(marker)].data(),
            scale,
            residual,
            static_cast<std::size_t>(fixture.samples));
    }
};

struct BedAvx2Ops
{
    const CompactGenotypeFixture& fixture;

    [[nodiscard]] auto dot(Eigen::Index marker, const double* residual) const
        -> double
    {
        return dot_bed_avx2(
            fixture.bed_column(marker),
            fixture.luts[static_cast<std::size_t>(marker)].data(),
            residual,
            fixture.samples);
    }

    auto axpy(Eigen::Index marker, double scale, double* residual) const -> void
    {
        axpy_bed_avx2(
            fixture.bed_column(marker),
            fixture.luts[static_cast<std::size_t>(marker)].data(),
            scale,
            residual,
            fixture.samples);
    }
};

#endif

template <typename Ops>
[[nodiscard]] auto run_dot_sweep(
    const CompactGenotypeFixture& fixture,
    const Ops& ops,
    const Eigen::VectorXd& residual) -> double
{
    double checksum = 0.0;
    for (Eigen::Index marker = 0; marker < fixture.markers; ++marker)
    {
        checksum += ops.dot(marker, residual.data());
    }
    return checksum;
}

template <typename Ops>
[[nodiscard]] auto run_update_sweep(
    const CompactGenotypeFixture& fixture,
    const Ops& ops,
    Eigen::VectorXd& residual) -> double
{
    double checksum = 0.0;
    for (Eigen::Index marker = 0; marker < fixture.markers; ++marker)
    {
        const double rhs = ops.dot(marker, residual.data());
        const double scale
            = fixture.alpha[static_cast<std::size_t>(marker)] + (1e-10 * rhs);
        ops.axpy(marker, scale, residual.data());
        checksum += rhs;
    }
    return checksum;
}

template <typename MaterializeDot>
[[nodiscard]] auto run_scratch_update_sweep(
    const CompactGenotypeFixture& fixture,
    MaterializeDot materialize_dot,
    Eigen::VectorXd& residual,
    Eigen::VectorXd& scratch) -> double
{
    double checksum = 0.0;
    for (Eigen::Index marker = 0; marker < fixture.markers; ++marker)
    {
        const double rhs
            = materialize_dot(marker, residual.data(), scratch.data());
        const double scale
            = fixture.alpha[static_cast<std::size_t>(marker)] + (1e-10 * rhs);
        residual.noalias() += scale * scratch;
        checksum += rhs;
    }
    return checksum;
}

template <typename Ops>
auto benchmark_dot(
    ankerl::nanobench::Bench& bench,
    std::string_view name,
    const CompactGenotypeFixture& fixture,
    const Ops& ops) -> void
{
    const std::string benchmark_name{name};
    bench.run(
        benchmark_name,
        [&]
        {
            const double checksum
                = run_dot_sweep(fixture, ops, fixture.initial_residual);
            ankerl::nanobench::doNotOptimizeAway(checksum);
        });
}

template <typename Ops>
auto benchmark_update(
    ankerl::nanobench::Bench& bench,
    std::string_view name,
    const CompactGenotypeFixture& fixture,
    const Ops& ops,
    Eigen::VectorXd& residual) -> void
{
    const std::string benchmark_name{name};
    bench.run(
        benchmark_name,
        [&]
        {
            residual = fixture.initial_residual;
            const double checksum = run_update_sweep(fixture, ops, residual);
            ankerl::nanobench::doNotOptimizeAway(checksum);
            ankerl::nanobench::doNotOptimizeAway(residual.data());
        });
}

template <typename MaterializeDot>
auto benchmark_scratch_update(
    ankerl::nanobench::Bench& bench,
    std::string_view name,
    const CompactGenotypeFixture& fixture,
    MaterializeDot materialize_dot,
    Eigen::VectorXd& residual,
    Eigen::VectorXd& scratch) -> void
{
    const std::string benchmark_name{name};
    bench.run(
        benchmark_name,
        [&]
        {
            residual = fixture.initial_residual;
            const double checksum = run_scratch_update_sweep(
                fixture, materialize_dot, residual, scratch);
            ankerl::nanobench::doNotOptimizeAway(checksum);
            ankerl::nanobench::doNotOptimizeAway(residual.data());
        });
}

auto check_scalar_kernels(const CompactGenotypeFixture& fixture) -> void
{
    Eigen::VectorXd expected(fixture.samples);
    Eigen::VectorXd actual(fixture.samples);
    Eigen::VectorXd scratch(fixture.samples);
    for (Eigen::Index marker = 0; marker < fixture.markers; ++marker)
    {
        const auto* lut = fixture.luts[static_cast<std::size_t>(marker)].data();
        const double dense_dot
            = fixture.dense.col(marker).dot(fixture.initial_residual);
        REQUIRE(
            gelex::bayes::detail::dot_scalar(
                fixture.uint8_column(marker),
                lut,
                fixture.initial_residual.data(),
                static_cast<std::size_t>(fixture.samples))
            == Catch::Approx(dense_dot).epsilon(1e-12));
        REQUIRE(
            dot_bed_scalar(
                fixture.bed_column(marker),
                lut,
                fixture.initial_residual.data(),
                fixture.samples)
            == Catch::Approx(dense_dot).epsilon(1e-12));

        expected = fixture.initial_residual;
        expected.noalias() += 0.125 * fixture.dense.col(marker);
        actual = fixture.initial_residual;
        gelex::bayes::detail::axpy_scalar(
            fixture.uint8_column(marker),
            lut,
            0.125,
            actual.data(),
            static_cast<std::size_t>(fixture.samples));
        REQUIRE(actual.isApprox(expected, 1e-14));
        actual = fixture.initial_residual;
        axpy_bed_scalar(
            fixture.bed_column(marker),
            lut,
            0.125,
            actual.data(),
            fixture.samples);
        REQUIRE(actual.isApprox(expected, 1e-14));

        const double uint8_scratch_dot = materialize_dot_uint8_scalar(
            fixture.uint8_column(marker),
            lut,
            fixture.initial_residual.data(),
            scratch.data(),
            fixture.samples);
        REQUIRE(scratch.isApprox(fixture.dense.col(marker), 1e-14));
        REQUIRE(uint8_scratch_dot == Catch::Approx(dense_dot).epsilon(1e-12));
        const double bed_scratch_dot = materialize_dot_bed_scalar(
            fixture.bed_column(marker),
            lut,
            fixture.initial_residual.data(),
            scratch.data(),
            fixture.samples);
        REQUIRE(scratch.isApprox(fixture.dense.col(marker), 1e-14));
        REQUIRE(bed_scratch_dot == Catch::Approx(dense_dot).epsilon(1e-12));
    }
}

#if GELEX_BENCHMARK_X86_64 && (defined(__GNUC__) || defined(__clang__))

auto check_avx2_kernels(const CompactGenotypeFixture& fixture) -> void
{
    Eigen::VectorXd expected(fixture.samples);
    Eigen::VectorXd actual(fixture.samples);
    Eigen::VectorXd scratch(fixture.samples);
    for (Eigen::Index marker = 0; marker < fixture.markers; ++marker)
    {
        const auto* lut = fixture.luts[static_cast<std::size_t>(marker)].data();
        const double dense_dot
            = fixture.dense.col(marker).dot(fixture.initial_residual);
        REQUIRE(
            gelex::bayes::detail::dot_avx2(
                fixture.uint8_column(marker),
                lut,
                fixture.initial_residual.data(),
                static_cast<std::size_t>(fixture.samples))
            == Catch::Approx(dense_dot).epsilon(1e-12));
        REQUIRE(
            dot_bed_avx2(
                fixture.bed_column(marker),
                lut,
                fixture.initial_residual.data(),
                fixture.samples)
            == Catch::Approx(dense_dot).epsilon(1e-12));

        expected = fixture.initial_residual;
        expected.noalias() += 0.125 * fixture.dense.col(marker);
        actual = fixture.initial_residual;
        gelex::bayes::detail::axpy_avx2(
            fixture.uint8_column(marker),
            lut,
            0.125,
            actual.data(),
            static_cast<std::size_t>(fixture.samples));
        REQUIRE(actual.isApprox(expected, 1e-14));
        actual = fixture.initial_residual;
        axpy_bed_avx2(
            fixture.bed_column(marker),
            lut,
            0.125,
            actual.data(),
            fixture.samples);
        REQUIRE(actual.isApprox(expected, 1e-14));

        const double uint8_scratch_dot = materialize_dot_uint8_avx2(
            fixture.uint8_column(marker),
            lut,
            fixture.initial_residual.data(),
            scratch.data(),
            fixture.samples);
        REQUIRE(scratch.isApprox(fixture.dense.col(marker), 1e-14));
        REQUIRE(uint8_scratch_dot == Catch::Approx(dense_dot).epsilon(1e-12));
        const double bed_scratch_dot = materialize_dot_bed_avx2(
            fixture.bed_column(marker),
            lut,
            fixture.initial_residual.data(),
            scratch.data(),
            fixture.samples);
        REQUIRE(scratch.isApprox(fixture.dense.col(marker), 1e-14));
        REQUIRE(bed_scratch_dot == Catch::Approx(dense_dot).epsilon(1e-12));
    }
}

auto check_avx512_kernels(const CompactGenotypeFixture& fixture) -> void
{
    Eigen::VectorXd expected(fixture.samples);
    Eigen::VectorXd actual(fixture.samples);
    for (Eigen::Index marker = 0; marker < fixture.markers; ++marker)
    {
        const auto* lut = fixture.luts[static_cast<std::size_t>(marker)].data();
        const double dense_dot
            = fixture.dense.col(marker).dot(fixture.initial_residual);
        REQUIRE(
            gelex::bayes::detail::dot_avx512(
                fixture.uint8_column(marker),
                lut,
                fixture.initial_residual.data(),
                static_cast<std::size_t>(fixture.samples))
            == Catch::Approx(dense_dot).epsilon(1e-12));

        expected = fixture.initial_residual;
        expected.noalias() += 0.125 * fixture.dense.col(marker);
        actual = fixture.initial_residual;
        gelex::bayes::detail::axpy_avx512(
            fixture.uint8_column(marker),
            lut,
            0.125,
            actual.data(),
            static_cast<std::size_t>(fixture.samples));
        REQUIRE(actual.isApprox(expected, 1e-14));
    }
}

#endif

}  // namespace

TEST_CASE(
    "compact genotype kernels agree with dense genotype",
    "[!benchmark][mcmc][compact_genotype][correctness]")
{
    const CompactGenotypeFixture fixture{517, 13};
    check_scalar_kernels(fixture);
#if GELEX_BENCHMARK_X86_64 && (defined(__GNUC__) || defined(__clang__))
    if (gelex::bayes::detail::supports_avx2())
    {
        check_avx2_kernels(fixture);
    }
    if (gelex::bayes::detail::supports_avx512())
    {
        check_avx512_kernels(fixture);
    }
#endif
}

TEST_CASE(
    "compact genotype dot sweep",
    "[!benchmark][mcmc][compact_genotype][dot]")
{
    constexpr std::array<Shape, 7> SHAPES{{
        {.samples = 512, .markers = 1024},
        {.samples = 512, .markers = 2048},
        {.samples = 512, .markers = 4096},
        {.samples = 512, .markers = 8192},
        {.samples = 512, .markers = 16384},
        {.samples = 2048, .markers = 4096},
        {.samples = 8192, .markers = 1024},
    }};

    for (const auto [samples, markers] : SHAPES)
    {
        const CompactGenotypeFixture fixture{samples, markers};
        const double entries
            = static_cast<double>(samples) * static_cast<double>(markers);
        ankerl::nanobench::Bench bench;
        bench
            .title(
                "compact genotype dot: n=" + std::to_string(samples)
                + ", p=" + std::to_string(markers))
            .unit("genotype")
            .batch(entries)
            .relative(true)
            .epochs(9)
            .warmup(1)
            .minEpochIterations(1)
            .minEpochTime(std::chrono::milliseconds{50});

        benchmark_dot(bench, "dense/current", fixture, DenseOps{fixture});
        benchmark_dot(bench, "uint8/scalar", fixture, Uint8ScalarOps{fixture});
        benchmark_dot(bench, "bed2/scalar", fixture, BedScalarOps{fixture});
#if GELEX_BENCHMARK_X86_64 && (defined(__GNUC__) || defined(__clang__))
        if (gelex::bayes::detail::supports_avx2())
        {
            benchmark_dot(
                bench, "uint8/avx2-permute", fixture, Uint8Avx2Ops{fixture});
            benchmark_dot(
                bench, "bed2/avx2-table-permute", fixture, BedAvx2Ops{fixture});
        }
        if (gelex::bayes::detail::supports_avx512())
        {
            benchmark_dot(
                bench,
                "uint8/avx512-permute",
                fixture,
                Uint8Avx512Ops{fixture});
        }
#endif
    }
}

TEST_CASE(
    "compact genotype dot and axpy sweep",
    "[!benchmark][mcmc][compact_genotype][update]")
{
    constexpr std::array<Shape, 10> SHAPES{{
        {.samples = 512, .markers = 1024},
        {.samples = 512, .markers = 2048},
        {.samples = 512, .markers = 4096},
        {.samples = 512, .markers = 8192},
        {.samples = 512, .markers = 16384},
        {.samples = 2048, .markers = 4096},
        {.samples = 8192, .markers = 1024},
        {.samples = 32768, .markers = 256},
        {.samples = 65536, .markers = 128},
        {.samples = 131072, .markers = 64},
    }};

    for (const auto [samples, markers] : SHAPES)
    {
        const CompactGenotypeFixture fixture{samples, markers};
        Eigen::VectorXd residual(samples);
        Eigen::VectorXd scratch(samples);
        const double entries
            = static_cast<double>(samples) * static_cast<double>(markers);
        ankerl::nanobench::Bench bench;
        bench
            .title(
                "compact genotype dot+axpy: n=" + std::to_string(samples)
                + ", p=" + std::to_string(markers))
            .unit("genotype")
            .batch(entries)
            .relative(true)
            .epochs(9)
            .warmup(1)
            .minEpochIterations(1)
            .minEpochTime(std::chrono::milliseconds{50});

        benchmark_update(
            bench, "dense/current", fixture, DenseOps{fixture}, residual);
        benchmark_update(
            bench,
            "uint8/scalar-direct",
            fixture,
            Uint8ScalarOps{fixture},
            residual);
        benchmark_update(
            bench,
            "bed2/scalar-direct",
            fixture,
            BedScalarOps{fixture},
            residual);
        benchmark_scratch_update(
            bench,
            "uint8/scalar-scratch",
            fixture,
            [&](Eigen::Index marker, const double* y, double* x)
            {
                return materialize_dot_uint8_scalar(
                    fixture.uint8_column(marker),
                    fixture.luts[static_cast<std::size_t>(marker)].data(),
                    y,
                    x,
                    fixture.samples);
            },
            residual,
            scratch);
        benchmark_scratch_update(
            bench,
            "bed2/scalar-scratch",
            fixture,
            [&](Eigen::Index marker, const double* y, double* x)
            {
                return materialize_dot_bed_scalar(
                    fixture.bed_column(marker),
                    fixture.luts[static_cast<std::size_t>(marker)].data(),
                    y,
                    x,
                    fixture.samples);
            },
            residual,
            scratch);

#if GELEX_BENCHMARK_X86_64 && (defined(__GNUC__) || defined(__clang__))
        if (gelex::bayes::detail::supports_avx2())
        {
            benchmark_update(
                bench,
                "uint8/avx2-direct",
                fixture,
                Uint8Avx2Ops{fixture},
                residual);
            benchmark_update(
                bench,
                "bed2/avx2-direct",
                fixture,
                BedAvx2Ops{fixture},
                residual);
            benchmark_scratch_update(
                bench,
                "uint8/avx2-scratch",
                fixture,
                [&](Eigen::Index marker, const double* y, double* x)
                {
                    return materialize_dot_uint8_avx2(
                        fixture.uint8_column(marker),
                        fixture.luts[static_cast<std::size_t>(marker)].data(),
                        y,
                        x,
                        fixture.samples);
                },
                residual,
                scratch);
            benchmark_scratch_update(
                bench,
                "bed2/avx2-scratch",
                fixture,
                [&](Eigen::Index marker, const double* y, double* x)
                {
                    return materialize_dot_bed_avx2(
                        fixture.bed_column(marker),
                        fixture.luts[static_cast<std::size_t>(marker)].data(),
                        y,
                        x,
                        fixture.samples);
                },
                residual,
                scratch);
        }
        if (gelex::bayes::detail::supports_avx512())
        {
            benchmark_update(
                bench,
                "uint8/avx512-direct",
                fixture,
                Uint8Avx512Ops{fixture},
                residual);
        }
#endif
    }
}

TEST_CASE(
    "register-direct compact kernels with 100k markers",
    "[!benchmark][mcmc][compact_genotype][direct][p100k]")
{
    constexpr std::array<Shape, 3> SHAPES{{
        {.samples = 512, .markers = 100000},
        {.samples = 2048, .markers = 100000},
        {.samples = 4096, .markers = 100000},
    }};

    if (!gelex::bayes::detail::supports_avx2())
    {
        SKIP("AVX2 and FMA are required");
    }
#if GELEX_BENCHMARK_X86_64 && (defined(__GNUC__) || defined(__clang__))
    for (const auto [samples, markers] : SHAPES)
    {
        const CompactGenotypeFixture fixture{samples, markers};
        const double entries
            = static_cast<double>(samples) * static_cast<double>(markers);

        ankerl::nanobench::Bench dot_bench;
        dot_bench
            .title(
                "register-direct dot: n=" + std::to_string(samples)
                + ", p=" + std::to_string(markers))
            .unit("genotype")
            .batch(entries)
            .relative(true)
            .epochs(9)
            .warmup(1)
            .minEpochIterations(1)
            .minEpochTime(std::chrono::milliseconds{50});
        benchmark_dot(dot_bench, "dense/current", fixture, DenseOps{fixture});
        benchmark_dot(dot_bench, "dense/cblas", fixture, BlasDenseOps{fixture});
        benchmark_dot(
            dot_bench,
            "uint8/avx2-register-direct",
            fixture,
            Uint8Avx2Ops{fixture});
        benchmark_dot(
            dot_bench,
            "bed2/avx2-register-direct",
            fixture,
            BedAvx2Ops{fixture});

        Eigen::VectorXd residual(samples);
        ankerl::nanobench::Bench update_bench;
        update_bench
            .title(
                "register-direct dot+axpy: n=" + std::to_string(samples)
                + ", p=" + std::to_string(markers))
            .unit("genotype")
            .batch(entries)
            .relative(true)
            .epochs(9)
            .warmup(1)
            .minEpochIterations(1)
            .minEpochTime(std::chrono::milliseconds{50});
        benchmark_update(
            update_bench,
            "dense/current",
            fixture,
            DenseOps{fixture},
            residual);
        benchmark_update(
            update_bench,
            "dense/cblas",
            fixture,
            BlasDenseOps{fixture},
            residual);
        benchmark_update(
            update_bench,
            "uint8/avx2-register-direct",
            fixture,
            Uint8Avx2Ops{fixture},
            residual);
        benchmark_update(
            update_bench,
            "bed2/avx2-register-direct",
            fixture,
            BedAvx2Ops{fixture},
            residual);
    }
#endif
}
