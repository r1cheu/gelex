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

#ifndef GELEX_BAYES_DETAIL_COMMON_KERNEL_H_
#define GELEX_BAYES_DETAIL_COMMON_KERNEL_H_

#include <Eigen/Core>
#include <random>

#include "gelex/bayes/design.h"
#include "gelex/bayes/state.h"
#include "gelex/bayes/stats/normal_sampler.h"
#include "gelex/bayes/stats/scaled_inv_chi2_sampler.h"
#include "gelex/bayes/variance/parameter.h"
#include "gelex/data/fixed_design.h"

namespace gelex::detail
{

class FixedEffectKernel
{
   public:
    auto step(
        const FixedDesign& design,
        FixedEffectState& state,
        ResidualState& residual,
        std::mt19937_64& rng) -> void
    {
        normal_.reset();
        for (Eigen::Index index = 0; index < state.coefficients.size(); ++index)
        {
            const auto column = design.X().col(index);
            const double old_value = state.coefficients(index);
            const double quadratic = design.xtx_diag()(index);
            const double linear = column.dot(residual.adjusted_response)
                                  + (quadratic * old_value);
            const double new_value = normal_.draw(
                {.mean = linear / quadratic,
                 .var = residual.variance / quadratic},
                rng);
            state.coefficients(index) = new_value;
            residual.adjusted_response.array()
                += (old_value - new_value) * column.array();
        }
    }

   private:
    NormalSampler<double> normal_{0.0};
};

class RandomEffectKernel
{
   public:
    explicit RandomEffectKernel(const VarianceParameter& parameter)
        : variance_sampler_{parameter.prior()}
    {
    }

    auto step(
        const bayes::RandomDesign& design,
        RandomEffectState& state,
        ResidualState& residual,
        std::mt19937_64& rng) -> void
    {
        normal_.reset();
        normal_.set_prior_var(state.variance);
        variance_sampler_.reset();

        state.fitted_values.setZero();
        for (Eigen::Index index = 0; index < state.coefficients.size(); ++index)
        {
            const auto column = design.X().col(index);
            const double old_value = state.coefficients(index);
            const double quadratic = design.xtx_diag()(index);
            const double linear = column.dot(residual.adjusted_response)
                                  + (quadratic * old_value);
            const double new_value = normal_(
                {.quadratic = quadratic,
                 .linear = linear,
                 .scale = residual.variance},
                rng);
            state.coefficients(index) = new_value;
            state.fitted_values.array() += new_value * column.array();
            const double difference = old_value - new_value;
            if (difference != 0.0)
            {
                residual.adjusted_response.array()
                    += difference * column.array();
            }
        }

        state.variance = variance_sampler_(
            {.n = state.coefficients.size(),
             .sum_squares = state.coefficients.squaredNorm()},
            rng);
    }

   private:
    ScaledInvChi2Sampler<double> variance_sampler_;
    NormalSampler<double> normal_{0.0};
};

class ResidualVarianceKernel
{
   public:
    explicit ResidualVarianceKernel(const VarianceParameter& parameter)
        : sampler_{parameter.prior()}
    {
    }

    auto step(ResidualState& state, std::mt19937_64& rng) -> void
    {
        sampler_.reset();
        state.variance = sampler_(
            {.n = state.adjusted_response.size(),
             .sum_squares = state.adjusted_response.squaredNorm()},
            rng);
    }

   private:
    ScaledInvChi2Sampler<double> sampler_;
};

}  // namespace gelex::detail

#endif  // GELEX_BAYES_DETAIL_COMMON_KERNEL_H_
