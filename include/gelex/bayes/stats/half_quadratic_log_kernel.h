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

#ifndef GELEX_BAYES_STATS_HALF_QUADRATIC_LOG_KERNEL_H_
#define GELEX_BAYES_STATS_HALF_QUADRATIC_LOG_KERNEL_H_

#include <array>
#include <cassert>
#include <cmath>
#include <numbers>
#include <random>

#include "gelex/bayes/stats/quadratic_log_kernel.h"
#include "gelex/bayes/stats/truncated_normal_distribution.h"
#include "gelex/infra/normal.h"

namespace gelex
{

class HalfQuadraticLogKernel
{
   public:
    class Evaluation
    {
       public:
        [[nodiscard]] auto log_integral(HalfLine support) const -> double
        {
            return support == HalfLine::Positive ? positive_log_integral_
                                                 : negative_log_integral_;
        }

        [[nodiscard]] auto truncated_normal_parameters(HalfLine support) const
            -> TruncatedNormalDistribution<>::param_type
        {
            return TruncatedNormalDistribution<>::param_type{
                normal_parameters_.mean(),
                normal_parameters_.stddev(),
                support};
        }

       private:
        friend class HalfQuadraticLogKernel;

        Evaluation(
            std::normal_distribution<double>::param_type normal_parameters,
            std::array<double, 2> log_integrals)
            : normal_parameters_{normal_parameters},
              negative_log_integral_{log_integrals[0]},
              positive_log_integral_{log_integrals[1]}
        {
        }

        std::normal_distribution<double>::param_type normal_parameters_;
        double negative_log_integral_;
        double positive_log_integral_;
    };

    explicit HalfQuadraticLogKernel(QuadraticLogKernel kernel) : kernel_{kernel}
    {
    }

    [[nodiscard]] auto evaluate() const -> Evaluation
    {
        const auto normal_parameters = kernel_.normal_parameters();
        const double standardized_mean
            = normal_parameters.mean() / normal_parameters.stddev();
        const double base_log_integral = kernel_.log_integral();
        return Evaluation{
            normal_parameters,
            std::array{
                base_log_integral + log_norm_cdf(-standardized_mean),
                base_log_integral + log_norm_cdf(standardized_mean)}};
    }

    [[nodiscard]] auto operator+(const QuadraticLogKernel& other) const
        -> HalfQuadraticLogKernel
    {
        return HalfQuadraticLogKernel{kernel_ + other};
    }

   private:
    QuadraticLogKernel kernel_;
};

[[nodiscard]] inline auto operator+(
    const QuadraticLogKernel& lhs,
    const HalfQuadraticLogKernel& rhs) -> HalfQuadraticLogKernel
{
    return rhs + lhs;
}

[[nodiscard]] inline auto make_half_normal_prior(double variance)
    -> HalfQuadraticLogKernel
{
    assert(variance > 0.0);

    const double precision = 1.0 / variance;
    const double constant
        = std::log(2.0)
          - (0.5 * (std::log(2.0 * std::numbers::pi) + std::log(variance)));
    return HalfQuadraticLogKernel{QuadraticLogKernel{precision, 0.0, constant}};
}

}  // namespace gelex

#endif  // GELEX_BAYES_STATS_HALF_QUADRATIC_LOG_KERNEL_H_
