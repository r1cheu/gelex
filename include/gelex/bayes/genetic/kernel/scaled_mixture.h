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

#ifndef GELEX_BAYES_GENETIC_KERNEL_SCALED_MIXTURE_H_
#define GELEX_BAYES_GENETIC_KERNEL_SCALED_MIXTURE_H_

#include <Eigen/Core>
#include <array>
#include <cstddef>
#include <cstdint>
#include <random>
#include <span>

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

template <std::size_t ClassCount, MixtureWeightUpdate WeightUpdate>
class ScaledMixtureKernel
{
    using Prior = ScaledMixturePrior<ClassCount, WeightUpdate>;
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
          allocation_updater_{prior.probabilities},
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
        variance_sampler_.reset();
        allocation_updater_.begin_sweep(family_state.probabilities);

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
                rng);
            const double new_value
                = sample.class_index == 0
                      ? 0.0
                      : normal_.draw(sample.coefficient_posterior, rng);

            coefficients(marker) = new_value;
            family_state.assignment(marker)
                = static_cast<std::uint8_t>(sample.class_index);

            std::array<bayes::AxpyTarget, 3> fitted_targets{};
            std::size_t fitted_target_count = 0;
            const double coefficient_delta = new_value - old_value;
            if (coefficient_delta != 0.0)
            {
                // NOLINTNEXTLINE(cppcoreguidelines-pro-bounds-constant-array-index)
                fitted_targets[fitted_target_count++] = bayes::AxpyTarget{
                    -coefficient_delta, residual.adjusted_response};
            }

            const auto append_component_target
                = [&](std::size_t class_index, double delta)
            {
                if (class_index == 0 || delta == 0.0)
                {
                    return;
                }

                const auto component_index
                    = static_cast<Eigen::Index>(class_index - 1);
                // NOLINTNEXTLINE(cppcoreguidelines-pro-bounds-constant-array-index)
                fitted_targets[fitted_target_count++] = bayes::AxpyTarget{
                    delta, family_state.fitted_values.col(component_index)};
            };

            if (old_class_index == sample.class_index)
            {
                append_component_target(old_class_index, coefficient_delta);
            }
            else
            {
                append_component_target(old_class_index, -old_value);
                append_component_target(sample.class_index, new_value);
            }

            if (fitted_target_count != 0)
            {
                projection.axpy(
                    marker,
                    std::span<const bayes::AxpyTarget>{
                        fitted_targets.data(), fitted_target_count});
            }

            if (sample.class_index != 0)
            {
                ++active_count;
                // NOLINTNEXTLINE(cppcoreguidelines-pro-bounds-constant-array-index)
                scaled_sum_squares
                    += (new_value * new_value) / scales_[sample.class_index];
            }
        }

        family_state.variance = variance_sampler_(
            {.n = active_count, .sum_squares = scaled_sum_squares}, rng);
        allocation_updater_.update(family_state.probabilities, rng);
    }

   private:
    auto draw_component(
        const NormalSampler<double>::Kernel& likelihood,
        double variance,
        std::mt19937_64& rng) -> ComponentSample
    {
        std::array<Posterior, ClassCount> posteriors{};
        std::array<double, ClassCount> log_likelihood_kernels{};
        for (std::size_t class_index = 1; class_index < ClassCount;
             ++class_index)
        {
            // NOLINTNEXTLINE(cppcoreguidelines-pro-bounds-constant-array-index)
            auto& coefficient_posterior = posteriors[class_index];
            // NOLINTNEXTLINE(cppcoreguidelines-pro-bounds-constant-array-index)
            const double class_scale = scales_[class_index];
            coefficient_posterior
                = normal_.set_prior_var(variance * class_scale)
                      .posterior_with_logL(likelihood);
            // NOLINTNEXTLINE(cppcoreguidelines-pro-bounds-constant-array-index)
            log_likelihood_kernels[class_index]
                = coefficient_posterior.log_likelihood_kernel;
        }

        const auto allocation_posterior
            = allocation_updater_.posterior(log_likelihood_kernels);
        const std::size_t class_index
            = allocation_updater_.draw(allocation_posterior, rng);
        // NOLINTNEXTLINE(cppcoreguidelines-pro-bounds-constant-array-index)
        return {
            .class_index = class_index,
            .coefficient_posterior = posteriors[class_index].params};
    }

    ScaledInvChi2Sampler<double> variance_sampler_;
    CategoricalAllocationUpdater<WeightUpdate, ClassCount> allocation_updater_;
    std::array<double, ClassCount> scales_;
    NormalSampler<double> normal_{0.0};
};

}  // namespace gelex::detail

#endif  // GELEX_BAYES_GENETIC_KERNEL_SCALED_MIXTURE_H_
