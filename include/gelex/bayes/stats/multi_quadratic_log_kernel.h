// Copyright 2026 RuLei Chen
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//     http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.

#ifndef GELEX_BAYES_STATS_MULTI_QUADRATIC_LOG_KERNEL_H_
#define GELEX_BAYES_STATS_MULTI_QUADRATIC_LOG_KERNEL_H_

#include <Eigen/Cholesky>
#include <Eigen/Core>
#include <cassert>
#include <cmath>
#include <numbers>
#include <utility>

namespace gelex
{

class MultiQuadraticLogKernel
{
   public:
    MultiQuadraticLogKernel(
        Eigen::Matrix2d quadratic,
        Eigen::Vector2d linear,
        double constant)
        : quadratic_{std::move(quadratic)},
          linear_{std::move(linear)},
          constant_{constant}
    {
        assert(quadratic_.allFinite());
        assert(linear_.allFinite());
        assert(std::isfinite(constant_));
        assert(quadratic_.isApprox(quadratic_.transpose()));
        assert(quadratic_.ldlt().isPositive());
    }

    [[nodiscard]] auto quadratic() const noexcept -> const Eigen::Matrix2d&
    {
        return quadratic_;
    }

    [[nodiscard]] auto linear() const noexcept -> const Eigen::Vector2d&
    {
        return linear_;
    }

    [[nodiscard]] auto constant() const noexcept -> double { return constant_; }

    [[nodiscard]] auto operator+(const MultiQuadraticLogKernel& other) const
        -> MultiQuadraticLogKernel
    {
        return MultiQuadraticLogKernel{
            quadratic_ + other.quadratic_,
            linear_ + other.linear_,
            constant_ + other.constant_};
    }

   private:
    Eigen::Matrix2d quadratic_;
    Eigen::Vector2d linear_;
    double constant_;
};

[[nodiscard]] inline auto make_multi_normal_prior(Eigen::Matrix2d covariance)
    -> MultiQuadraticLogKernel
{
    assert(covariance.allFinite());
    assert(covariance.isApprox(covariance.transpose()));

    const Eigen::LLT<Eigen::Matrix2d> covariance_factor{covariance};
    assert(covariance_factor.info() == Eigen::Success);
    const Eigen::Matrix2d precision
        = covariance_factor.solve(Eigen::Matrix2d::Identity());
    const double log_determinant
        = 2.0 * covariance_factor.matrixLLT().diagonal().array().log().sum();
    const double constant
        = -std::log(2.0 * std::numbers::pi) - (0.5 * log_determinant);
    return MultiQuadraticLogKernel{
        precision, Eigen::Vector2d::Zero(), constant};
}

}  // namespace gelex

#endif  // GELEX_BAYES_STATS_MULTI_QUADRATIC_LOG_KERNEL_H_
