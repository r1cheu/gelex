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

#ifndef GELEX_BAYES_GENETIC_KERNEL_SPIKE_SLAB_H_
#define GELEX_BAYES_GENETIC_KERNEL_SPIKE_SLAB_H_

#include <Eigen/Core>
#include <cstdint>
#include <random>

#include "gelex/bayes/genetic/kernel/allocation_updater.h"
#include "gelex/bayes/genetic/prior.h"
#include "gelex/bayes/genetic/state.h"
#include "gelex/bayes/genotype/design.h"
#include "gelex/bayes/state.h"
#include "gelex/bayes/stats/normal_sampler.h"
#include "gelex/bayes/stats/scaled_inv_chi2_sampler.h"
#include "gelex/genetic_mode.h"

namespace gelex::detail
{

template <VarianceLayout Kind, MixtureWeightUpdate WeightUpdate>
class SpikeSlabKernel
{
    using Prior = SpikeSlabPrior<Kind, WeightUpdate>;
    using State = GeneticModeState<SpikeSlabState<Kind>>;

   public:
    explicit SpikeSlabKernel(const Prior& prior)
        : variance_sampler_{prior.variance.prior()},
          allocation_updater_{prior.probability}
    {
    }

    template <GeneticMode Mode>
    auto step(
        const bayes::GeneticDesign& design,
        State& state,
        ResidualState& residual,
        std::mt19937_64& rng) -> void
    {
        const auto& projection = design.projection(Mode);
        const auto valid_indices = projection.valid_indices();
        const auto& xtx_diag = projection.xtx_diag();
        auto& coefficients = state.coefficients;
        auto& family_state = state.family_state;
        auto& variance = family_state.variance;

        previous_adjusted_response_ = residual.adjusted_response;
        normal_.reset();
        variance_sampler_.reset();
        allocation_updater_.begin_sweep(family_state.probability);
        if constexpr (Kind == VarianceLayout::Pooled)
        {
            normal_.set_prior_var(variance);
        }
        Eigen::Index pooled_active_count = 0;
        double pooled_sum_squares = 0.0;
        for (const Eigen::Index marker : valid_indices)
        {
            if constexpr (Kind == VarianceLayout::Unpooled)
            {
                normal_.set_prior_var(variance(marker));
            }

            const double old_value = coefficients(marker);
            const double linear
                = projection.dot(marker, residual.adjusted_response)
                  + (xtx_diag(marker) * old_value);
            const auto coefficient_posterior = normal_.posterior_with_logL(
                NormalSampler<double>::Kernel{
                    .quadratic = xtx_diag(marker),
                    .linear = linear,
                    .scale = residual.variance});
            const auto allocation_posterior = allocation_updater_.posterior(
                {0.0, coefficient_posterior.log_likelihood_kernel});
            const bool is_active
                = allocation_updater_.draw(allocation_posterior, rng);
            const double new_value
                = is_active ? normal_.draw(coefficient_posterior.params, rng)
                            : 0.0;
            const double squared_effect = new_value * new_value;

            coefficients(marker) = new_value;
            family_state.assignment(marker)
                = static_cast<std::uint8_t>(is_active);
            projection.axpy(
                marker, old_value - new_value, residual.adjusted_response);

            if constexpr (Kind == VarianceLayout::Pooled)
            {
                if (is_active)
                {
                    ++pooled_active_count;
                    pooled_sum_squares += squared_effect;
                }
            }
            if constexpr (Kind == VarianceLayout::Unpooled)
            {
                // The inactive slab coefficient is integrated out, so it is
                // not a zero-valued Gaussian observation of its variance.
                variance(marker) = variance_sampler_(
                    {.n = is_active ? 1 : 0, .sum_squares = squared_effect},
                    rng);
            }
        }

        state.family_state.fitted_values.col(0).noalias()
            += previous_adjusted_response_ - residual.adjusted_response;
        if constexpr (Kind == VarianceLayout::Pooled)
        {
            variance = variance_sampler_(
                {.n = pooled_active_count, .sum_squares = pooled_sum_squares},
                rng);
        }
        allocation_updater_.update(family_state.probability, rng);
    }

   private:
    ScaledInvChi2Sampler<double> variance_sampler_;
    BinaryAllocationUpdater<WeightUpdate> allocation_updater_;
    NormalSampler<double> normal_{0.0};
    Eigen::VectorXd previous_adjusted_response_;
};

}  // namespace gelex::detail

#endif  // GELEX_BAYES_GENETIC_KERNEL_SPIKE_SLAB_H_
