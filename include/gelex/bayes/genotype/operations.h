// Copyright 2026 RuLei Chen
// SPDX-License-Identifier: Apache-2.0

#ifndef GELEX_BAYES_GENOTYPE_OPERATIONS_H_
#define GELEX_BAYES_GENOTYPE_OPERATIONS_H_

#include <Eigen/Core>
#include <cstdint>
#include <span>

namespace gelex::bayes
{

struct AxpyTarget
{
    AxpyTarget() = default;

    AxpyTarget(double scale, Eigen::Ref<Eigen::VectorXd> target) noexcept
        : scale{scale}, target{target}
    {
    }

    double scale{};
    std::span<double> target;
};

[[nodiscard]] auto dot(
    std::span<const std::uint8_t> genotype_column,
    const Eigen::Ref<const Eigen::Array4d>& lut,
    std::span<const double> rhs) noexcept -> double;

auto multiply(
    std::span<const std::uint8_t> genotype_column,
    const Eigen::Ref<const Eigen::Array4d>& lut,
    double scale,
    std::span<double> target) noexcept -> void;

auto axpy(
    std::span<const std::uint8_t> genotype_column,
    const Eigen::Ref<const Eigen::Array4d>& lut,
    double scale,
    std::span<double> target) noexcept -> void;

auto axpy(
    std::span<const std::uint8_t> genotype_column,
    const Eigen::Ref<const Eigen::Array4d>& lut,
    std::span<const AxpyTarget> targets) noexcept -> void;

}  // namespace gelex::bayes

#endif  // GELEX_BAYES_GENOTYPE_OPERATIONS_H_
