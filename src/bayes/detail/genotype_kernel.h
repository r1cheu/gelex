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

#ifndef GELEX_BAYES_DETAIL_GENOTYPE_KERNEL_H_
#define GELEX_BAYES_DETAIL_GENOTYPE_KERNEL_H_

#include <cstddef>
#include <cstdint>
#include <span>

#include "gelex/bayes/genotype_kernel.h"

namespace gelex::bayes::detail
{

using DotKernel = double (*)(
    const std::uint8_t* genotype_column,
    const double* lut,
    const double* residual,
    std::size_t size) noexcept;

using AxpyKernel = void (*)(
    const std::uint8_t* genotype_column,
    const double* lut,
    double scale,
    double* residual,
    std::size_t size) noexcept;

using MultiTargetAxpyKernel = void (*)(
    const std::uint8_t* genotype_column,
    const double* lut,
    std::span<const AxpyTarget> targets,
    std::size_t size) noexcept;

[[nodiscard]] auto dot_scalar(
    const std::uint8_t* genotype_column,
    const double* lut,
    const double* residual,
    std::size_t size) noexcept -> double;

auto axpy_scalar(
    const std::uint8_t* genotype_column,
    const double* lut,
    double scale,
    double* residual,
    std::size_t size) noexcept -> void;

auto axpy_multi_target_scalar(
    const std::uint8_t* genotype_column,
    const double* lut,
    std::span<const AxpyTarget> targets,
    std::size_t size) noexcept -> void;

[[nodiscard]] auto dot_avx2(
    const std::uint8_t* genotype_column,
    const double* lut,
    const double* residual,
    std::size_t size) noexcept -> double;

auto axpy_avx2(
    const std::uint8_t* genotype_column,
    const double* lut,
    double scale,
    double* residual,
    std::size_t size) noexcept -> void;

auto axpy_multi_target_avx2(
    const std::uint8_t* genotype_column,
    const double* lut,
    std::span<const AxpyTarget> targets,
    std::size_t size) noexcept -> void;

[[nodiscard]] auto dot_avx512(
    const std::uint8_t* genotype_column,
    const double* lut,
    const double* residual,
    std::size_t size) noexcept -> double;

auto axpy_avx512(
    const std::uint8_t* genotype_column,
    const double* lut,
    double scale,
    double* residual,
    std::size_t size) noexcept -> void;

auto axpy_multi_target_avx512(
    const std::uint8_t* genotype_column,
    const double* lut,
    std::span<const AxpyTarget> targets,
    std::size_t size) noexcept -> void;

[[nodiscard]] auto supports_avx2() noexcept -> bool;
[[nodiscard]] auto supports_avx512() noexcept -> bool;
[[nodiscard]] auto select_dot_kernel() noexcept -> DotKernel;
[[nodiscard]] auto select_axpy_kernel() noexcept -> AxpyKernel;
[[nodiscard]] auto select_multi_target_axpy_kernel() noexcept
    -> MultiTargetAxpyKernel;

}  // namespace gelex::bayes::detail

#endif  // GELEX_BAYES_DETAIL_GENOTYPE_KERNEL_H_
