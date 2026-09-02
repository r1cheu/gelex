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

#ifndef GELEX_BAYES_GENETIC_SCALED_MIXTURE_KERNEL_H_
#define GELEX_BAYES_GENETIC_SCALED_MIXTURE_KERNEL_H_

#include <Eigen/Core>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <random>

#include "gelex/bayes/design_data.h"
#include "gelex/bayes/detail/fitted_value_update.h"
#include "gelex/bayes/genetic/prior.h"
#include "gelex/bayes/genetic/probability_updater.h"
#include "gelex/bayes/genetic/state.h"
#include "gelex/bayes/state.h"
#include "gelex/infra/stats/normal_sampler.h"
#include "gelex/infra/stats/scaled_inv_chi2_sampler.h"
#include "gelex/types/genetic_mode.h"

namespace gelex::detail
{

template <std::size_t ClassCount, UpdatePolicy ProbabilitiesUpdate>
class ScaledMixtureKernel
{
    using Prior = ScaledMixturePrior<ClassCount, ProbabilitiesUpdate>;
    using State = GeneticModeState<ScaledMixtureState<ClassCount>>;
    using Posterior = NormalSampler<double>::Posterior;
    using PosteriorParams = NormalSampler<double>::Params;

    struct ComponentSample
    {
        std::size_t class_index{};
        PosteriorParams coefficient_posterior{};
    };

   public:
    explicit ScaledMixtureKernel(const Prior& prior)
        : variance_sampler_{prior.variance.prior()},
          probabilities_updater_{prior.probabilities},
          scales_{prior.scales}
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

        normal_.reset();
        uniform_.reset();
        variance_sampler_.reset();
        if constexpr (ProbabilitiesUpdate == UpdatePolicy::Sampled)
        {
            probabilities_updater_.reset();
        }

        std::array<double, ClassCount> log_probabilities{};
        for (std::size_t class_index = 0; class_index < ClassCount;
             ++class_index)
        {
            log_probabilities[class_index]
                = std::log(family_state.probabilities[class_index]);
        }

        Eigen::Index active_count = 0;
        double scaled_sum_squares = 0.0;
        for (const Eigen::Index marker : valid_indices)
        {
            const double old_value = coefficients(marker);
            const auto old_class_index
                = static_cast<std::size_t>(family_state.assignment(marker));
            const double linear
                = projection.dot(marker, residual.adjusted_response)
                  + (xtx_diag(marker) * old_value);
            const auto sample = draw_component(
                {.quadratic = xtx_diag(marker),
                 .linear = linear,
                 .scale = residual.variance},
                family_state.variance,
                log_probabilities,
                rng);
            const double new_value
                = sample.class_index == 0
                      ? 0.0
                      : normal_.draw(sample.coefficient_posterior, rng);

            coefficients(marker) = new_value;
            family_state.assignment(marker)
                = static_cast<std::uint8_t>(sample.class_index);

            apply_fitted_value_transition(
                projection,
                marker,
                FittedValueTransition{
                    .old_class_index = old_class_index,
                    .new_class_index = sample.class_index,
                    .old_coefficient = old_value,
                    .new_coefficient = new_value},
                family_state,
                residual.adjusted_response);

            if (sample.class_index != 0)
            {
                ++active_count;
                scaled_sum_squares
                    += (new_value * new_value) / scales_[sample.class_index];
            }
            if constexpr (ProbabilitiesUpdate == UpdatePolicy::Sampled)
            {
                probabilities_updater_.observe(sample.class_index);
            }
        }

        family_state.variance = variance_sampler_(
            {.n = active_count, .sum_squares = scaled_sum_squares}, rng);
        if constexpr (ProbabilitiesUpdate == UpdatePolicy::Sampled)
        {
            probabilities_updater_.update(family_state.probabilities, rng);
        }
    }

   private:
    auto draw_component(
        const NormalSampler<double>::Kernel& likelihood,
        double variance,
        const std::array<double, ClassCount>& log_probabilities,
        std::mt19937_64& rng) -> ComponentSample
    {
        std::array<Posterior, ClassCount> posteriors{};
        std::array<double, ClassCount> weights{};
        weights[0] = log_probabilities[0];
        double max_log_weight = weights[0];
        for (std::size_t class_index = 1; class_index < ClassCount;
             ++class_index)
        {
            posteriors[class_index]
                = normal_.set_prior_var(variance * scales_[class_index])
                      .posterior_with_logL(likelihood);
            weights[class_index] = posteriors[class_index].log_likelihood_kernel
                                   + log_probabilities[class_index];
            max_log_weight = std::max(max_log_weight, weights[class_index]);
        }

        double total_weight = 0.0;
        for (double& weight : weights)
        {
            weight = std::exp(weight - max_log_weight);
            total_weight += weight;
        }

        const double threshold = uniform_(rng) * total_weight;
        double cumulative_weight = 0.0;
        for (std::size_t class_index = 0; class_index < ClassCount;
             ++class_index)
        {
            cumulative_weight += weights[class_index];
            if (threshold < cumulative_weight)
            {
                return {
                    .class_index = class_index,
                    .coefficient_posterior = posteriors[class_index].params};
            }
        }
        return {
            .class_index = ClassCount - 1,
            .coefficient_posterior = posteriors.back().params};
    }

    ScaledInvChi2Sampler<double> variance_sampler_;
    [[no_unique_address]]
    SimplexUpdater<ProbabilitiesUpdate, ClassCount> probabilities_updater_;
    std::array<double, ClassCount> scales_;
    NormalSampler<double> normal_{0.0};
    std::uniform_real_distribution<double> uniform_{0.0, 1.0};
};

}  // namespace gelex::detail

#endif  // GELEX_BAYES_GENETIC_SCALED_MIXTURE_KERNEL_H_
