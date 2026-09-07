// Copyright 2026 RuLei Chen
// SPDX-License-Identifier: Apache-2.0

#ifndef GELEX_BAYES_GENOTYPE_DETAIL_OPERATIONS_H_
#define GELEX_BAYES_GENOTYPE_DETAIL_OPERATIONS_H_

#include <cstddef>
#include <cstdint>
#include <span>

#include "gelex/bayes/genotype/operations.h"

namespace gelex::bayes::detail
{

using DotImpl = double (*)(
    const std::uint8_t* genotype_column,
    const double* lut,
    const double* rhs,
    std::size_t size) noexcept;

using MultiplyImpl = void (*)(
    const std::uint8_t* genotype_column,
    const double* lut,
    double scale,
    double* target,
    std::size_t size) noexcept;

using AxpyImpl = void (*)(
    const std::uint8_t* genotype_column,
    const double* lut,
    double scale,
    double* target,
    std::size_t size) noexcept;

using MultiTargetAxpyImpl = void (*)(
    const std::uint8_t* genotype_column,
    const double* lut,
    std::span<const AxpyTarget> targets,
    std::size_t size) noexcept;

[[nodiscard]] auto dot_scalar(
    const std::uint8_t* genotype_column,
    const double* lut,
    const double* rhs,
    std::size_t size) noexcept -> double;

auto multiply_scalar(
    const std::uint8_t* genotype_column,
    const double* lut,
    double scale,
    double* target,
    std::size_t size) noexcept -> void;

auto axpy_scalar(
    const std::uint8_t* genotype_column,
    const double* lut,
    double scale,
    double* target,
    std::size_t size) noexcept -> void;

auto axpy_multi_target_scalar(
    const std::uint8_t* genotype_column,
    const double* lut,
    std::span<const AxpyTarget> targets,
    std::size_t size) noexcept -> void;

[[nodiscard]] auto dot_avx2(
    const std::uint8_t* genotype_column,
    const double* lut,
    const double* rhs,
    std::size_t size) noexcept -> double;

auto multiply_avx2(
    const std::uint8_t* genotype_column,
    const double* lut,
    double scale,
    double* target,
    std::size_t size) noexcept -> void;

auto axpy_avx2(
    const std::uint8_t* genotype_column,
    const double* lut,
    double scale,
    double* target,
    std::size_t size) noexcept -> void;

auto axpy_multi_target_avx2(
    const std::uint8_t* genotype_column,
    const double* lut,
    std::span<const AxpyTarget> targets,
    std::size_t size) noexcept -> void;

[[nodiscard]] auto dot_avx512(
    const std::uint8_t* genotype_column,
    const double* lut,
    const double* rhs,
    std::size_t size) noexcept -> double;

auto multiply_avx512(
    const std::uint8_t* genotype_column,
    const double* lut,
    double scale,
    double* target,
    std::size_t size) noexcept -> void;

auto axpy_avx512(
    const std::uint8_t* genotype_column,
    const double* lut,
    double scale,
    double* target,
    std::size_t size) noexcept -> void;

auto axpy_multi_target_avx512(
    const std::uint8_t* genotype_column,
    const double* lut,
    std::span<const AxpyTarget> targets,
    std::size_t size) noexcept -> void;

[[nodiscard]] auto supports_avx2() noexcept -> bool;
[[nodiscard]] auto supports_avx512() noexcept -> bool;
[[nodiscard]] auto select_dot_impl() noexcept -> DotImpl;
[[nodiscard]] auto select_multiply_impl() noexcept -> MultiplyImpl;
[[nodiscard]] auto select_axpy_impl() noexcept -> AxpyImpl;
[[nodiscard]] auto select_multi_target_axpy_impl() noexcept
    -> MultiTargetAxpyImpl;

}  // namespace gelex::bayes::detail

#endif  // GELEX_BAYES_GENOTYPE_DETAIL_OPERATIONS_H_
