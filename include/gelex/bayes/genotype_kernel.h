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

#ifndef GELEX_BAYES_GENOTYPE_KERNEL_H_
#define GELEX_BAYES_GENOTYPE_KERNEL_H_

#include <Eigen/Core>
#include <cstddef>
#include <cstdint>
#include <span>

namespace gelex::bayes
{

struct AxpyTarget
{
    AxpyTarget() = default;

    AxpyTarget(double scale, Eigen::Ref<Eigen::VectorXd> values) noexcept
        : scale{scale},
          values{values.data(), static_cast<std::size_t>(values.size())}
    {
    }

    double scale{};
    std::span<double> values;
};

[[nodiscard]] auto dot(
    std::span<const std::uint8_t> genotype_column,
    const Eigen::Ref<const Eigen::Array4d>& lut,
    std::span<const double> residual) noexcept -> double;

auto axpy(
    std::span<const std::uint8_t> genotype_column,
    const Eigen::Ref<const Eigen::Array4d>& lut,
    double scale,
    std::span<double> residual) noexcept -> void;

auto axpy(
    std::span<const std::uint8_t> genotype_column,
    const Eigen::Ref<const Eigen::Array4d>& lut,
    std::span<const AxpyTarget> targets) noexcept -> void;

}  // namespace gelex::bayes

#endif  // GELEX_BAYES_GENOTYPE_KERNEL_H_
