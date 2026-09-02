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
#include <cstddef>
#include <random>

#include "gelex/bayes/design.h"
#include "gelex/bayes/detail/normal_variance_conjugate_updater.h"
#include "gelex/bayes/parameter.h"
#include "gelex/bayes/state.h"
#include "gelex/bayes/stats/quadratic_log_kernel.h"
#include "gelex/data/fixed_design.h"

namespace gelex::detail
{

inline auto update_fixed_effects(
    const FixedDesign& design,
    FixedEffectState& state,
    ResidualState& residual,
    std::mt19937_64& rng) -> void
{
    std::normal_distribution<double> normal_distribution;
    for (Eigen::Index index = 0; index < state.coefficients.size(); ++index)
    {
        const auto column = design.X().col(index);
        const double old_value = state.coefficients(index);
        const double quadratic = design.xtx_diag()(index);
        const double linear
            = column.dot(residual.adjusted_response) + (quadratic * old_value);
        const auto likelihood = gelex::make_coefficient_likelihood(
            quadratic, linear, residual.variance);
        const double new_value
            = normal_distribution(rng, likelihood.normal_parameters());
        state.coefficients(index) = new_value;
        residual.adjusted_response.array()
            += (old_value - new_value) * column.array();
    }
}

class RandomEffectKernel
{
   public:
    explicit RandomEffectKernel(const VarianceParameter& parameter)
        : variance_updater_{parameter.prior}
    {
    }

    auto step(
        const bayes::RandomDesign& design,
        RandomEffectState& state,
        ResidualState& residual,
        std::mt19937_64& rng) -> void
    {
        std::normal_distribution<double> normal_distribution;
        const auto coefficient_prior = gelex::make_normal_prior(state.variance);

        state.fitted_values.setZero();
        for (Eigen::Index index = 0; index < state.coefficients.size(); ++index)
        {
            const auto column = design.X().col(index);
            const double old_value = state.coefficients(index);
            const double quadratic = design.xtx_diag()(index);
            const double linear = column.dot(residual.adjusted_response)
                                  + (quadratic * old_value);
            const auto likelihood = gelex::make_coefficient_likelihood(
                quadratic, linear, residual.variance);
            const auto posterior = likelihood + coefficient_prior;
            const double new_value
                = normal_distribution(rng, posterior.normal_parameters());
            state.coefficients(index) = new_value;
            state.fitted_values.array() += new_value * column.array();
            const double difference = old_value - new_value;
            if (difference != 0.0)
            {
                residual.adjusted_response.array()
                    += difference * column.array();
            }
        }

        variance_updater_.update(
            state.variance,
            static_cast<std::size_t>(state.coefficients.size()),
            state.coefficients.squaredNorm(),
            rng);
    }

   private:
    NormalVarianceConjugateUpdater variance_updater_;
};

class ResidualVarianceKernel
{
   public:
    explicit ResidualVarianceKernel(const VarianceParameter& parameter)
        : variance_updater_{parameter.prior}
    {
    }

    auto step(ResidualState& state, std::mt19937_64& rng) -> void
    {
        variance_updater_.update(
            state.variance,
            static_cast<std::size_t>(state.adjusted_response.size()),
            state.adjusted_response.squaredNorm(),
            rng);
    }

   private:
    NormalVarianceConjugateUpdater variance_updater_;
};

}  // namespace gelex::detail

#endif  // GELEX_BAYES_DETAIL_COMMON_KERNEL_H_
