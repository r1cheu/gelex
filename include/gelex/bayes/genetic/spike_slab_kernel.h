// Copyright 2026 RuLei Chen
// SPDX-License-Identifier: Apache-2.0

#ifndef GELEX_BAYES_GENETIC_SPIKE_SLAB_KERNEL_H_
#define GELEX_BAYES_GENETIC_SPIKE_SLAB_KERNEL_H_

#include <Eigen/Core>
#include <array>
#include <cstddef>
#include <cstdint>
#include <random>
#include <type_traits>
#include <variant>

#include "gelex/bayes/detail/normal_variance_conjugate_updater.h"
#include "gelex/bayes/genetic/detail/coefficient_likelihood.h"
#include "gelex/bayes/genetic/detail/dirichlet_conjugate_updater.h"
#include "gelex/bayes/genetic/detail/normal_prior_provider.h"
#include "gelex/bayes/genetic/spike_slab.h"
#include "gelex/bayes/genetic/types.h"
#include "gelex/bayes/genotype/design.h"
#include "gelex/bayes/state.h"
#include "gelex/bayes/stats/log_categorical_distribution.h"
#include "gelex/genetic_mode.h"

namespace gelex
{

template <VarianceLayout Kind, MixtureWeightUpdate WeightUpdate>
class SpikeSlabKernel
{
    using prior_type = SpikeSlabPrior<Kind, WeightUpdate>;
    using state_type = SpikeSlabState<Kind, WeightUpdate>;
    using probability_updater_type = std::conditional_t<
        WeightUpdate == MixtureWeightUpdate::Enabled,
        detail::DirichletConjugateUpdater<2>,
        std::monostate>;

   public:
    explicit SpikeSlabKernel(const prior_type& prior)
        : variance_updater_{prior.variance.prior},
          probability_updater_{
              [&]() -> probability_updater_type
              {
                  if constexpr (WeightUpdate == MixtureWeightUpdate::Enabled)
                  {
                      return probability_updater_type{prior.probability.prior};
                  }
                  else
                  {
                      return {};
                  }
              }()}
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
        auto& variance = state.variance();

        previous_adjusted_response_ = residual.adjusted_response;
        std::normal_distribution<double> normal_distribution;
        const auto log_probabilities = make_log_weights(
            std::array{1.0 - state.probability(), state.probability()});
        const auto normal_prior_for_marker
            = detail::make_normal_prior_provider<Kind>(variance);

        double pooled_sum_squares = 0.0;
        for (const Eigen::Index marker : valid_indices)
        {
            const double old_value = coefficients(marker);
            const auto likelihood = detail::make_coefficient_likelihood(
                projection, marker, old_value, residual);
            const auto slab_kernel
                = likelihood + normal_prior_for_marker(marker);

            const auto allocation_parameters = make_mixture_posterior_weights(
                log_probabilities, std::array{0.0, slab_kernel.log_integral()});
            const std::size_t allocation
                = allocation_distribution_(rng, allocation_parameters);
            const bool is_active = allocation == 1;
            const double new_value
                = is_active ? normal_distribution(
                                  rng, slab_kernel.normal_parameters())
                            : 0.0;
            const double squared_effect = new_value * new_value;

            state.transition(marker, new_value, is_active);
            projection.axpy(
                marker, old_value - new_value, residual.adjusted_response);

            if constexpr (Kind == VarianceLayout::Pooled)
            {
                if (is_active)
                {
                    pooled_sum_squares += squared_effect;
                }
            }
            if constexpr (Kind == VarianceLayout::Unpooled)
            {
                // An inactive slab is integrated out and contributes no
                // Gaussian observation to its variance posterior.
                variance_updater_.update(
                    variance(marker), is_active ? 1 : 0, squared_effect, rng);
            }
        }

        previous_adjusted_response_ -= residual.adjusted_response;
        state.transition(previous_adjusted_response_);
        if constexpr (Kind == VarianceLayout::Pooled)
        {
            variance_updater_.update(
                variance, state.class_counts()[1], pooled_sum_squares, rng);
        }
        if constexpr (WeightUpdate == MixtureWeightUpdate::Enabled)
        {
            auto allocation_counts = state.class_counts();
            // Skipped markers remain NULL but do not contribute to the
            // posterior.
            allocation_counts[0]
                -= static_cast<std::size_t>(coefficients.size())
                   - valid_indices.size();
            state.set_probability(
                probability_updater_.draw(allocation_counts, rng)[1]);
        }
    }

   private:
    detail::NormalVarianceConjugateUpdater variance_updater_;
    [[no_unique_address]] probability_updater_type probability_updater_;
    LogCategoricalDistribution<2> allocation_distribution_;
    Eigen::VectorXd previous_adjusted_response_;
};

template <VarianceLayout Kind, MixtureWeightUpdate WeightUpdate>
[[nodiscard]] auto make_kernel(const SpikeSlabPrior<Kind, WeightUpdate>& prior)
{
    return SpikeSlabKernel<Kind, WeightUpdate>{prior};
}

}  // namespace gelex

#endif  // GELEX_BAYES_GENETIC_SPIKE_SLAB_KERNEL_H_
