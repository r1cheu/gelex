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

#include "gelex/bayes/detail/normal_variance_conjugate_updater.h"
#include "gelex/bayes/genetic/detail/coefficient_likelihood.h"
#include "gelex/bayes/genetic/detail/dirichlet_conjugate_updater.h"
#include "gelex/bayes/genetic/prior.h"
#include "gelex/bayes/genetic/state.h"
#include "gelex/bayes/genotype/design.h"
#include "gelex/bayes/state.h"
#include "gelex/bayes/stats/log_categorical_distribution.h"
#include "gelex/bayes/stats/quadratic_log_kernel.h"
#include "gelex/genetic_mode.h"

namespace gelex::detail
{

template <std::size_t ClassCount, MixtureWeightUpdate WeightUpdate>
class ScaledMixtureKernel
{
    using Prior = ScaledMixturePrior<ClassCount, WeightUpdate>;
    using State = GeneticModeState<ScaledMixtureState<ClassCount>>;
    using CoefficientParameters = std::normal_distribution<double>::param_type;

    struct ComponentSample
    {
        std::size_t class_index{};
        CoefficientParameters coefficient_parameters;
    };

   public:
    explicit ScaledMixtureKernel(const Prior& prior)
        : variance_updater_{prior.variance.prior},
          probability_updater_{make_dirichlet_conjugate_updater<ClassCount>(
              prior.probabilities)},
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
        auto& coefficients = state.coefficients;
        auto& family_state = state.family_state;

        std::normal_distribution<double> normal_distribution;
        const auto log_probabilities
            = make_log_weights(family_state.probabilities);
        std::array<std::size_t, ClassCount> allocation_counts{};

        std::size_t active_count = 0;
        double scaled_sum_squares = 0.0;
        for (const Eigen::Index marker : valid_indices)
        {
            const double old_value = coefficients(marker);
            const auto old_class_index
                = static_cast<std::size_t>(family_state.assignment(marker));
            const auto likelihood = make_coefficient_likelihood(
                projection, marker, old_value, residual);
            const auto sample = draw_component(
                likelihood, family_state.variance, log_probabilities, rng);
            if constexpr (WeightUpdate == MixtureWeightUpdate::Enabled)
            {
                // NOLINTNEXTLINE(cppcoreguidelines-pro-bounds-constant-array-index)
                ++allocation_counts[sample.class_index];
            }
            const double new_value
                = sample.class_index == 0
                      ? 0.0
                      : normal_distribution(rng, sample.coefficient_parameters);

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

        variance_updater_.update(
            family_state.variance, active_count, scaled_sum_squares, rng);
        probability_updater_.update(
            family_state.probabilities, allocation_counts, rng);
    }

   private:
    auto draw_component(
        const QuadraticLogKernel& likelihood,
        double variance,
        const std::array<double, ClassCount>& log_probabilities,
        std::mt19937_64& rng) -> ComponentSample
    {
        std::array<CoefficientParameters, ClassCount> coefficient_parameters{};
        std::array<double, ClassCount> component_log_integrals{};
        for (std::size_t class_index = 1; class_index < ClassCount;
             ++class_index)
        {
            // NOLINTNEXTLINE(cppcoreguidelines-pro-bounds-constant-array-index)
            const double class_scale = scales_[class_index];
            const auto coefficient_posterior
                = likelihood + make_normal_prior(variance * class_scale);
            // NOLINTNEXTLINE(cppcoreguidelines-pro-bounds-constant-array-index)
            coefficient_parameters[class_index]
                = coefficient_posterior.normal_parameters();
            // NOLINTNEXTLINE(cppcoreguidelines-pro-bounds-constant-array-index)
            component_log_integrals[class_index]
                = coefficient_posterior.log_integral();
        }

        const auto allocation_parameters = make_mixture_posterior_weights(
            log_probabilities, component_log_integrals);
        const std::size_t class_index
            = allocation_distribution_(rng, allocation_parameters);
        // NOLINTNEXTLINE(cppcoreguidelines-pro-bounds-constant-array-index)
        return {
            .class_index = class_index,
            .coefficient_parameters = coefficient_parameters[class_index]};
    }

    NormalVarianceConjugateUpdater variance_updater_;
    [[no_unique_address]] DirichletConjugateUpdater<ClassCount, WeightUpdate>
        probability_updater_;
    LogCategoricalDistribution<ClassCount> allocation_distribution_;
    std::array<double, ClassCount> scales_;
};

}  // namespace gelex::detail

#endif  // GELEX_BAYES_GENETIC_KERNEL_SCALED_MIXTURE_H_
