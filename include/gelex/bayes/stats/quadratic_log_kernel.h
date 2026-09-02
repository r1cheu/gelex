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

#ifndef GELEX_BAYES_STATS_QUADRATIC_LOG_KERNEL_H_
#define GELEX_BAYES_STATS_QUADRATIC_LOG_KERNEL_H_

#include <cassert>
#include <cmath>
#include <numbers>
#include <random>

namespace gelex
{

class QuadraticLogKernel
{
   public:
    // NOLINTNEXTLINE(bugprone-easily-swappable-parameters)
    QuadraticLogKernel(double quadratic, double linear, double constant)
        : quadratic_{quadratic}, linear_{linear}, constant_{constant}
    {
        assert(quadratic_ >= 0.0);
    }

    [[nodiscard]] auto normal_parameters() const
        -> std::normal_distribution<double>::param_type
    {
        assert(quadratic_ > 0.0);
        return std::normal_distribution<double>::param_type{
            linear_ / quadratic_, 1.0 / std::sqrt(quadratic_)};
    }

    [[nodiscard]] auto log_integral() const -> double
    {
        assert(quadratic_ > 0.0);
        const double mean = linear_ / quadratic_;
        return constant_ + (0.5 * linear_ * mean)
               + (0.5
                  * (std::log(2.0 * std::numbers::pi) - std::log(quadratic_)));
    }

    [[nodiscard]] auto operator+(const QuadraticLogKernel& other) const
        -> QuadraticLogKernel
    {
        return QuadraticLogKernel(
            quadratic_ + other.quadratic_,
            linear_ + other.linear_,
            constant_ + other.constant_);
    }

   private:
    double quadratic_;
    double linear_;
    double constant_;
};

// NOLINTBEGIN(bugprone-easily-swappable-parameters)
[[nodiscard]] inline auto make_coefficient_likelihood(
    double quadratic,
    double linear,
    double residual_variance) -> QuadraticLogKernel
{
    assert(quadratic >= 0.0);
    assert(residual_variance > 0.0);

    const double residual_precision = 1.0 / residual_variance;
    return QuadraticLogKernel{
        quadratic * residual_precision, linear * residual_precision, 0.0};
}
// NOLINTEND(bugprone-easily-swappable-parameters)

[[nodiscard]] inline auto make_normal_prior(double variance)
    -> QuadraticLogKernel
{
    assert(variance > 0.0);

    const double precision = 1.0 / variance;
    const double constant
        = -0.5 * (std::log(2.0 * std::numbers::pi) + std::log(variance));

    return QuadraticLogKernel{precision, 0.0, constant};
}

}  // namespace gelex

#endif  // GELEX_BAYES_STATS_QUADRATIC_LOG_KERNEL_H_
