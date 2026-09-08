// Copyright 2026 RuLei Chen
// SPDX-License-Identifier: Apache-2.0

#ifndef GELEX_BAYES_GENETIC_SCALED_MIXTURE_KERNEL_H_
#define GELEX_BAYES_GENETIC_SCALED_MIXTURE_KERNEL_H_

#include <Eigen/Core>
#include <array>
#include <cstddef>
#include <cstdint>
#include <random>
#include <span>
#include <type_traits>
#include <variant>

#include "gelex/bayes/detail/normal_variance_conjugate_updater.h"
#include "gelex/bayes/genetic/detail/apply_fitted_update.h"
#include "gelex/bayes/genetic/detail/coefficient_likelihood.h"
#include "gelex/bayes/genetic/detail/dirichlet_conjugate_updater.h"
#include "gelex/bayes/genetic/scaled_mixture.h"
#include "gelex/bayes/genetic/types.h"
#include "gelex/bayes/genotype/design.h"
#include "gelex/bayes/genotype/operations.h"
#include "gelex/bayes/state.h"
#include "gelex/bayes/stats/log_categorical_distribution.h"
#include "gelex/bayes/stats/quadratic_log_kernel.h"
#include "gelex/genetic_mode.h"

namespace gelex
{

template <MixtureWeightUpdate WeightUpdate>
class ScaledMixtureKernel
{
    using prior_type = ScaledMixturePrior<WeightUpdate>;
    using state_type = ScaledMixtureState<WeightUpdate>;
    using coefficient_parameters_type
        = std::normal_distribution<double>::param_type;

    static constexpr std::size_t class_count = state_type::class_count;
    using probability_updater_type = std::conditional_t<
        WeightUpdate == MixtureWeightUpdate::Enabled,
        detail::DirichletConjugateUpdater<class_count>,
        std::monostate>;

    struct ComponentSample
    {
        std::size_t class_index{};
        coefficient_parameters_type coefficient_parameters;
    };

   public:
    explicit ScaledMixtureKernel(const prior_type& prior)
        : variance_updater_{prior.variance.prior},
          probability_updater_{
              [&]() -> probability_updater_type
              {
                  if constexpr (WeightUpdate == MixtureWeightUpdate::Enabled)
                  {
                      return probability_updater_type{
                          prior.probabilities.prior};
                  }
                  else
                  {
                      return {};
                  }
              }()},
          scales_{prior.scales}
    {
    }

    template <GeneticMode Mode>
    auto step(
        const bayes::GeneticDesign& design,
        state_type& state,
        ResidualState& residual,
        std::mt19937_64& rng) -> void
    {
        const auto& projection = design.projection(Mode);
        const auto valid_indices = projection.valid_indices();
        const auto& coefficients = state.coefficients();

        std::normal_distribution<double> normal_distribution;
        const auto log_probabilities = make_log_weights(state.probabilities());

        double scaled_sum_squares = 0.0;
        for (const Eigen::Index marker : valid_indices)
        {
            const double old_value = coefficients(marker);
            const auto likelihood = detail::make_coefficient_likelihood(
                projection, marker, old_value, residual);
            const auto sample = draw_component(
                likelihood, state.variance(), log_probabilities, rng);
            const double new_value
                = sample.class_index == 0
                      ? 0.0
                      : normal_distribution(rng, sample.coefficient_parameters);

            const std::array extra_targets{bayes::AxpyTarget{
                old_value - new_value, residual.adjusted_response}};
            detail::apply_fitted_update(
                projection,
                marker,
                state.transition(
                    marker,
                    new_value,
                    static_cast<std::uint8_t>(sample.class_index)),
                std::span{extra_targets});

            if (sample.class_index != 0)
            {
                // NOLINTNEXTLINE(cppcoreguidelines-pro-bounds-constant-array-index)
                scaled_sum_squares
                    += (new_value * new_value) / scales_[sample.class_index];
            }
        }

        const auto active_count = static_cast<std::size_t>(coefficients.size())
                                  - state.class_counts()[0];
        variance_updater_.update(
            state.variance(), active_count, scaled_sum_squares, rng);
        if constexpr (WeightUpdate == MixtureWeightUpdate::Enabled)
        {
            auto allocation_counts = state.class_counts();
            // Skipped markers remain NULL but do not contribute to the
            // posterior.
            allocation_counts[0]
                -= static_cast<std::size_t>(coefficients.size())
                   - valid_indices.size();
            state.set_probabilities(
                probability_updater_.draw(allocation_counts, rng));
        }
    }

   private:
    auto draw_component(
        const QuadraticLogKernel& likelihood,
        double variance,
        const std::array<double, class_count>& log_probabilities,
        std::mt19937_64& rng) -> ComponentSample
    {
        std::array<coefficient_parameters_type, class_count>
            coefficient_parameters{};
        std::array<double, class_count> component_log_integrals{};
        for (std::size_t class_index = 1; class_index < class_count;
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

    detail::NormalVarianceConjugateUpdater variance_updater_;
    [[no_unique_address]] probability_updater_type probability_updater_;
    LogCategoricalDistribution<class_count> allocation_distribution_;
    std::array<double, class_count> scales_;
};

template <MixtureWeightUpdate WeightUpdate>
[[nodiscard]] auto make_kernel(const ScaledMixturePrior<WeightUpdate>& prior)
{
    return ScaledMixtureKernel<WeightUpdate>{prior};
}

}  // namespace gelex

#endif  // GELEX_BAYES_GENETIC_SCALED_MIXTURE_KERNEL_H_
