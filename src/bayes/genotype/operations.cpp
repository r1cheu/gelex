// Copyright 2026 RuLei Chen
// SPDX-License-Identifier: Apache-2.0

#include "gelex/bayes/genotype/operations.h"

#include <Eigen/Core>
#include <algorithm>
#include <cassert>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <span>

#include "bayes/genotype/detail/operations.h"

#if (defined(__x86_64__) || defined(_M_X64)) \
    && (defined(__GNUC__) || defined(__clang__))
#include <immintrin.h>
#define GELEX_X86_FUNCTION_TARGETS 1
#define GELEX_AVX2_TARGET __attribute__((target("avx2,fma")))
#define GELEX_AVX512_TARGET __attribute__((target("avx512f,avx512bw")))
#else
#define GELEX_X86_FUNCTION_TARGETS 0
#define GELEX_AVX2_TARGET
#define GELEX_AVX512_TARGET
#endif

namespace gelex::bayes::detail
{

auto dot_scalar(
    const std::uint8_t* genotype_column,
    const double* lut,
    const double* rhs,
    std::size_t size) noexcept -> double
{
    double sum = 0.0;
    for (std::size_t index = 0; index < size; ++index)
    {
        sum += lut[genotype_column[index]] * rhs[index];
    }
    return sum;
}

auto multiply_scalar(
    const std::uint8_t* genotype_column,
    const double* lut,
    double scale,
    double* target,
    std::size_t size) noexcept -> void
{
    for (std::size_t index = 0; index < size; ++index)
    {
        target[index] = scale * lut[genotype_column[index]];
    }
}

auto axpy_scalar(
    const std::uint8_t* genotype_column,
    const double* lut,
    double scale,
    double* target,
    std::size_t size) noexcept -> void
{
    for (std::size_t index = 0; index < size; ++index)
    {
        target[index] += scale * lut[genotype_column[index]];
    }
}

auto axpy_multi_target_scalar(
    const std::uint8_t* genotype_column,
    const double* lut,
    std::span<const AxpyTarget> targets,
    std::size_t size) noexcept -> void
{
    for (std::size_t index = 0; index < size; ++index)
    {
        const double value = lut[genotype_column[index]];
        for (const auto& update : targets)
        {
            update.target[index] += update.scale * value;
        }
    }
}

#if GELEX_X86_FUNCTION_TARGETS

GELEX_AVX2_TARGET static inline auto lookup_avx2(
    std::uint32_t four_codes,
    __m256i lut_parts) noexcept -> __m256d
{
    const __m128i code_bytes
        = _mm_cvtsi32_si128(static_cast<std::int32_t>(four_codes));
    const __m128i codes = _mm_cvtepu8_epi32(code_bytes);
    const __m128i even_indices = _mm_slli_epi32(codes, 1);
    const __m128i odd_indices = _mm_add_epi32(even_indices, _mm_set1_epi32(1));
    const __m256i indices = _mm256_set_m128i(
        _mm_unpackhi_epi32(even_indices, odd_indices),
        _mm_unpacklo_epi32(even_indices, odd_indices));
    return _mm256_castsi256_pd(_mm256_permutevar8x32_epi32(lut_parts, indices));
}

GELEX_AVX2_TARGET static inline auto load_values_avx2(
    const std::uint8_t* genotype_column,
    __m256i lut_parts) noexcept -> __m256d
{
    std::uint32_t four_codes = 0;
    std::memcpy(&four_codes, genotype_column, sizeof(four_codes));
    return lookup_avx2(four_codes, lut_parts);
}

GELEX_AVX2_TARGET static inline auto horizontal_sum_avx2(
    __m256d values) noexcept -> double
{
    const __m128d low = _mm256_castpd256_pd128(values);
    const __m128d high = _mm256_extractf128_pd(values, 1);
    const __m128d pair = _mm_add_pd(low, high);
    const __m128d swapped = _mm_shuffle_pd(pair, pair, 0x01);
    return _mm_cvtsd_f64(_mm_add_sd(pair, swapped));
}

GELEX_AVX2_TARGET auto dot_avx2(
    const std::uint8_t* genotype_column,
    const double* lut,
    const double* rhs,
    std::size_t size) noexcept -> double
{
    const __m256i lut_parts = _mm256_castpd_si256(_mm256_loadu_pd(lut));
    __m256d sum0 = _mm256_setzero_pd();
    __m256d sum1 = _mm256_setzero_pd();
    __m256d sum2 = _mm256_setzero_pd();
    __m256d sum3 = _mm256_setzero_pd();

    std::size_t index = 0;
    for (; index + 16 <= size; index += 16)
    {
        sum0 = _mm256_fmadd_pd(
            load_values_avx2(genotype_column + index, lut_parts),
            _mm256_loadu_pd(rhs + index),
            sum0);
        sum1 = _mm256_fmadd_pd(
            load_values_avx2(genotype_column + index + 4, lut_parts),
            _mm256_loadu_pd(rhs + index + 4),
            sum1);
        sum2 = _mm256_fmadd_pd(
            load_values_avx2(genotype_column + index + 8, lut_parts),
            _mm256_loadu_pd(rhs + index + 8),
            sum2);
        sum3 = _mm256_fmadd_pd(
            load_values_avx2(genotype_column + index + 12, lut_parts),
            _mm256_loadu_pd(rhs + index + 12),
            sum3);
    }

    __m256d sum
        = _mm256_add_pd(_mm256_add_pd(sum0, sum1), _mm256_add_pd(sum2, sum3));
    for (; index + 4 <= size; index += 4)
    {
        sum = _mm256_fmadd_pd(
            load_values_avx2(genotype_column + index, lut_parts),
            _mm256_loadu_pd(rhs + index),
            sum);
    }

    double result = horizontal_sum_avx2(sum);
    for (; index < size; ++index)
    {
        result += lut[genotype_column[index]] * rhs[index];
    }
    return result;
}

GELEX_AVX2_TARGET auto multiply_avx2(
    const std::uint8_t* genotype_column,
    const double* lut,
    double scale,
    double* target,
    std::size_t size) noexcept -> void
{
    const __m256i lut_parts = _mm256_castpd_si256(_mm256_loadu_pd(lut));
    const __m256d scale_values = _mm256_set1_pd(scale);

    std::size_t index = 0;
    for (; index + 4 <= size; index += 4)
    {
        _mm256_storeu_pd(
            target + index,
            _mm256_mul_pd(
                scale_values,
                load_values_avx2(genotype_column + index, lut_parts)));
    }
    for (; index < size; ++index)
    {
        target[index] = scale * lut[genotype_column[index]];
    }
}

GELEX_AVX2_TARGET auto axpy_avx2(
    const std::uint8_t* genotype_column,
    const double* lut,
    double scale,
    double* target,
    std::size_t size) noexcept -> void
{
    const __m256i lut_parts = _mm256_castpd_si256(_mm256_loadu_pd(lut));
    const __m256d scale_values = _mm256_set1_pd(scale);

    std::size_t index = 0;
    for (; index + 4 <= size; index += 4)
    {
        const __m256d updated = _mm256_fmadd_pd(
            scale_values,
            load_values_avx2(genotype_column + index, lut_parts),
            _mm256_loadu_pd(target + index));
        _mm256_storeu_pd(target + index, updated);
    }
    for (; index < size; ++index)
    {
        target[index] += scale * lut[genotype_column[index]];
    }
}

GELEX_AVX2_TARGET auto axpy_multi_target_avx2(
    const std::uint8_t* genotype_column,
    const double* lut,
    std::span<const AxpyTarget> targets,
    std::size_t size) noexcept -> void
{
    const __m256i lut_parts = _mm256_castpd_si256(_mm256_loadu_pd(lut));

    std::size_t index = 0;
    for (; index + 16 <= size; index += 16)
    {
        const __m256d values0
            = load_values_avx2(genotype_column + index, lut_parts);
        const __m256d values1
            = load_values_avx2(genotype_column + index + 4, lut_parts);
        const __m256d values2
            = load_values_avx2(genotype_column + index + 8, lut_parts);
        const __m256d values3
            = load_values_avx2(genotype_column + index + 12, lut_parts);
        for (const auto& update : targets)
        {
            const __m256d scale = _mm256_set1_pd(update.scale);
            _mm256_storeu_pd(
                update.target.data() + index,
                _mm256_fmadd_pd(
                    scale,
                    values0,
                    _mm256_loadu_pd(update.target.data() + index)));
            _mm256_storeu_pd(
                update.target.data() + index + 4,
                _mm256_fmadd_pd(
                    scale,
                    values1,
                    _mm256_loadu_pd(update.target.data() + index + 4)));
            _mm256_storeu_pd(
                update.target.data() + index + 8,
                _mm256_fmadd_pd(
                    scale,
                    values2,
                    _mm256_loadu_pd(update.target.data() + index + 8)));
            _mm256_storeu_pd(
                update.target.data() + index + 12,
                _mm256_fmadd_pd(
                    scale,
                    values3,
                    _mm256_loadu_pd(update.target.data() + index + 12)));
        }
    }
    for (; index + 4 <= size; index += 4)
    {
        const __m256d values
            = load_values_avx2(genotype_column + index, lut_parts);
        for (const auto& update : targets)
        {
            const __m256d updated = _mm256_fmadd_pd(
                _mm256_set1_pd(update.scale),
                values,
                _mm256_loadu_pd(update.target.data() + index));
            _mm256_storeu_pd(update.target.data() + index, updated);
        }
    }
    for (; index < size; ++index)
    {
        const double value = lut[genotype_column[index]];
        for (const auto& update : targets)
        {
            update.target[index] += update.scale * value;
        }
    }
}

GELEX_AVX512_TARGET static inline auto load_values_avx512(
    const std::uint8_t* genotype_column,
    __m512d lut_values) noexcept -> __m512d
{
    std::uint64_t eight_codes = 0;
    std::memcpy(&eight_codes, genotype_column, sizeof(eight_codes));
    const __m128i code_bytes
        = _mm_cvtsi64_si128(static_cast<std::int64_t>(eight_codes));
    const __m512i indices = _mm512_cvtepu8_epi64(code_bytes);
    return _mm512_permutexvar_pd(indices, lut_values);
}

GELEX_AVX512_TARGET auto dot_avx512(
    const std::uint8_t* genotype_column,
    const double* lut,
    const double* rhs,
    std::size_t size) noexcept -> double
{
    const __m512d lut_values = _mm512_setr_pd(
        lut[0], lut[1], lut[2], lut[3], lut[0], lut[1], lut[2], lut[3]);
    __m512d sum0 = _mm512_setzero_pd();
    __m512d sum1 = _mm512_setzero_pd();
    __m512d sum2 = _mm512_setzero_pd();
    __m512d sum3 = _mm512_setzero_pd();

    std::size_t index = 0;
    for (; index + 32 <= size; index += 32)
    {
        sum0 = _mm512_fmadd_pd(
            load_values_avx512(genotype_column + index, lut_values),
            _mm512_loadu_pd(rhs + index),
            sum0);
        sum1 = _mm512_fmadd_pd(
            load_values_avx512(genotype_column + index + 8, lut_values),
            _mm512_loadu_pd(rhs + index + 8),
            sum1);
        sum2 = _mm512_fmadd_pd(
            load_values_avx512(genotype_column + index + 16, lut_values),
            _mm512_loadu_pd(rhs + index + 16),
            sum2);
        sum3 = _mm512_fmadd_pd(
            load_values_avx512(genotype_column + index + 24, lut_values),
            _mm512_loadu_pd(rhs + index + 24),
            sum3);
    }

    __m512d sum
        = _mm512_add_pd(_mm512_add_pd(sum0, sum1), _mm512_add_pd(sum2, sum3));
    for (; index + 8 <= size; index += 8)
    {
        sum = _mm512_fmadd_pd(
            load_values_avx512(genotype_column + index, lut_values),
            _mm512_loadu_pd(rhs + index),
            sum);
    }

    double result = _mm512_reduce_add_pd(sum);
    for (; index < size; ++index)
    {
        result += lut[genotype_column[index]] * rhs[index];
    }
    return result;
}

GELEX_AVX512_TARGET auto multiply_avx512(
    const std::uint8_t* genotype_column,
    const double* lut,
    double scale,
    double* target,
    std::size_t size) noexcept -> void
{
    const __m512d lut_values = _mm512_setr_pd(
        lut[0], lut[1], lut[2], lut[3], lut[0], lut[1], lut[2], lut[3]);
    const __m512d scale_values = _mm512_set1_pd(scale);

    std::size_t index = 0;
    for (; index + 8 <= size; index += 8)
    {
        _mm512_storeu_pd(
            target + index,
            _mm512_mul_pd(
                scale_values,
                load_values_avx512(genotype_column + index, lut_values)));
    }
    for (; index < size; ++index)
    {
        target[index] = scale * lut[genotype_column[index]];
    }
}

GELEX_AVX512_TARGET auto axpy_avx512(
    const std::uint8_t* genotype_column,
    const double* lut,
    double scale,
    double* target,
    std::size_t size) noexcept -> void
{
    const __m512d lut_values = _mm512_setr_pd(
        lut[0], lut[1], lut[2], lut[3], lut[0], lut[1], lut[2], lut[3]);
    const __m512d scale_values = _mm512_set1_pd(scale);

    std::size_t index = 0;
    for (; index + 8 <= size; index += 8)
    {
        const __m512d updated = _mm512_fmadd_pd(
            scale_values,
            load_values_avx512(genotype_column + index, lut_values),
            _mm512_loadu_pd(target + index));
        _mm512_storeu_pd(target + index, updated);
    }
    for (; index < size; ++index)
    {
        target[index] += scale * lut[genotype_column[index]];
    }
}

GELEX_AVX512_TARGET auto axpy_multi_target_avx512(
    const std::uint8_t* genotype_column,
    const double* lut,
    std::span<const AxpyTarget> targets,
    std::size_t size) noexcept -> void
{
    const __m512d lut_values = _mm512_setr_pd(
        lut[0], lut[1], lut[2], lut[3], lut[0], lut[1], lut[2], lut[3]);

    std::size_t index = 0;
    for (; index + 32 <= size; index += 32)
    {
        const __m512d values0
            = load_values_avx512(genotype_column + index, lut_values);
        const __m512d values1
            = load_values_avx512(genotype_column + index + 8, lut_values);
        const __m512d values2
            = load_values_avx512(genotype_column + index + 16, lut_values);
        const __m512d values3
            = load_values_avx512(genotype_column + index + 24, lut_values);
        for (const auto& update : targets)
        {
            const __m512d scale = _mm512_set1_pd(update.scale);
            _mm512_storeu_pd(
                update.target.data() + index,
                _mm512_fmadd_pd(
                    scale,
                    values0,
                    _mm512_loadu_pd(update.target.data() + index)));
            _mm512_storeu_pd(
                update.target.data() + index + 8,
                _mm512_fmadd_pd(
                    scale,
                    values1,
                    _mm512_loadu_pd(update.target.data() + index + 8)));
            _mm512_storeu_pd(
                update.target.data() + index + 16,
                _mm512_fmadd_pd(
                    scale,
                    values2,
                    _mm512_loadu_pd(update.target.data() + index + 16)));
            _mm512_storeu_pd(
                update.target.data() + index + 24,
                _mm512_fmadd_pd(
                    scale,
                    values3,
                    _mm512_loadu_pd(update.target.data() + index + 24)));
        }
    }
    for (; index + 8 <= size; index += 8)
    {
        const __m512d values
            = load_values_avx512(genotype_column + index, lut_values);
        for (const auto& update : targets)
        {
            const __m512d updated = _mm512_fmadd_pd(
                _mm512_set1_pd(update.scale),
                values,
                _mm512_loadu_pd(update.target.data() + index));
            _mm512_storeu_pd(update.target.data() + index, updated);
        }
    }
    for (; index < size; ++index)
    {
        const double value = lut[genotype_column[index]];
        for (const auto& update : targets)
        {
            update.target[index] += update.scale * value;
        }
    }
}

auto supports_avx2() noexcept -> bool
{
    __builtin_cpu_init();
    return __builtin_cpu_supports("avx2") && __builtin_cpu_supports("fma");
}

auto supports_avx512() noexcept -> bool
{
    __builtin_cpu_init();
    return __builtin_cpu_supports("avx512f")
           && __builtin_cpu_supports("avx512bw");
}

#else

auto dot_avx2(
    const std::uint8_t* genotype_column,
    const double* lut,
    const double* rhs,
    std::size_t size) noexcept -> double
{
    return dot_scalar(genotype_column, lut, rhs, size);
}

auto multiply_avx2(
    const std::uint8_t* genotype_column,
    const double* lut,
    double scale,
    double* target,
    std::size_t size) noexcept -> void
{
    multiply_scalar(genotype_column, lut, scale, target, size);
}

auto axpy_avx2(
    const std::uint8_t* genotype_column,
    const double* lut,
    double scale,
    double* target,
    std::size_t size) noexcept -> void
{
    axpy_scalar(genotype_column, lut, scale, target, size);
}

auto axpy_multi_target_avx2(
    const std::uint8_t* genotype_column,
    const double* lut,
    std::span<const AxpyTarget> targets,
    std::size_t size) noexcept -> void
{
    axpy_multi_target_scalar(genotype_column, lut, targets, size);
}

auto dot_avx512(
    const std::uint8_t* genotype_column,
    const double* lut,
    const double* rhs,
    std::size_t size) noexcept -> double
{
    return dot_scalar(genotype_column, lut, rhs, size);
}

auto multiply_avx512(
    const std::uint8_t* genotype_column,
    const double* lut,
    double scale,
    double* target,
    std::size_t size) noexcept -> void
{
    multiply_scalar(genotype_column, lut, scale, target, size);
}

auto axpy_avx512(
    const std::uint8_t* genotype_column,
    const double* lut,
    double scale,
    double* target,
    std::size_t size) noexcept -> void
{
    axpy_scalar(genotype_column, lut, scale, target, size);
}

auto axpy_multi_target_avx512(
    const std::uint8_t* genotype_column,
    const double* lut,
    std::span<const AxpyTarget> targets,
    std::size_t size) noexcept -> void
{
    axpy_multi_target_scalar(genotype_column, lut, targets, size);
}

auto supports_avx2() noexcept -> bool
{
    return false;
}

auto supports_avx512() noexcept -> bool
{
    return false;
}

#endif

namespace
{

template <typename Impl>
[[nodiscard]] auto
select_impl(Impl scalar_impl, Impl avx2_impl, Impl avx512_impl) noexcept -> Impl
{
    if (supports_avx512())
    {
        return avx512_impl;
    }
    if (supports_avx2())
    {
        return avx2_impl;
    }
    return scalar_impl;
}

}  // namespace

auto select_dot_impl() noexcept -> DotImpl
{
    return select_impl<DotImpl>(dot_scalar, dot_avx2, dot_avx512);
}

auto select_multiply_impl() noexcept -> MultiplyImpl
{
    return select_impl<MultiplyImpl>(
        multiply_scalar, multiply_avx2, multiply_avx512);
}

auto select_axpy_impl() noexcept -> AxpyImpl
{
    return select_impl<AxpyImpl>(axpy_scalar, axpy_avx2, axpy_avx512);
}

auto select_multi_target_axpy_impl() noexcept -> MultiTargetAxpyImpl
{
    return select_impl<MultiTargetAxpyImpl>(
        axpy_multi_target_scalar,
        axpy_multi_target_avx2,
        axpy_multi_target_avx512);
}

}  // namespace gelex::bayes::detail

namespace gelex::bayes
{

auto dot(
    std::span<const std::uint8_t> genotype_column,
    const Eigen::Ref<const Eigen::Array4d>& lut,
    std::span<const double> rhs) noexcept -> double
{
    assert(genotype_column.size() == rhs.size());
    static const detail::DotImpl impl = detail::select_dot_impl();
    return impl(genotype_column.data(), lut.data(), rhs.data(), rhs.size());
}

auto multiply(
    std::span<const std::uint8_t> genotype_column,
    const Eigen::Ref<const Eigen::Array4d>& lut,
    double scale,
    std::span<double> target) noexcept -> void
{
    assert(genotype_column.size() == target.size());
    static const detail::MultiplyImpl impl = detail::select_multiply_impl();
    impl(
        genotype_column.data(),
        lut.data(),
        scale,
        target.data(),
        target.size());
}

auto axpy(
    std::span<const std::uint8_t> genotype_column,
    const Eigen::Ref<const Eigen::Array4d>& lut,
    double scale,
    std::span<double> target) noexcept -> void
{
    assert(genotype_column.size() == target.size());
    static const detail::AxpyImpl impl = detail::select_axpy_impl();
    impl(
        genotype_column.data(),
        lut.data(),
        scale,
        target.data(),
        target.size());
}

auto axpy(
    std::span<const std::uint8_t> genotype_column,
    const Eigen::Ref<const Eigen::Array4d>& lut,
    std::span<const AxpyTarget> targets) noexcept -> void
{
    assert(
        std::ranges::all_of(
            targets,
            [size = genotype_column.size()](const AxpyTarget& update)
            { return update.target.size() == size; }));
    static const detail::MultiTargetAxpyImpl impl
        = detail::select_multi_target_axpy_impl();
    impl(genotype_column.data(), lut.data(), targets, genotype_column.size());
}

}  // namespace gelex::bayes
