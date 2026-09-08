// Copyright 2026 RuLei Chen
// SPDX-License-Identifier: Apache-2.0

#ifndef GELEX_BAYES_GENETIC_GAUSSIAN_KERNEL_H_
#define GELEX_BAYES_GENETIC_GAUSSIAN_KERNEL_H_

#include <Eigen/Core>
#include <random>

#include "gelex/bayes/detail/normal_variance_conjugate_updater.h"
#include "gelex/bayes/genetic/detail/coefficient_likelihood.h"
#include "gelex/bayes/genetic/detail/normal_prior_provider.h"
#include "gelex/bayes/genetic/gaussian.h"
#include "gelex/bayes/genetic/types.h"
#include "gelex/bayes/genotype/design.h"
#include "gelex/bayes/genotype/projection.h"
#include "gelex/bayes/state.h"
#include "gelex/genetic_mode.h"

namespace gelex
{

template <VarianceLayout Kind>
class GaussianKernel
{
    using prior_type = GaussianPrior<Kind>;
    using state_type = GaussianState<Kind>;

   public:
    explicit GaussianKernel(const prior_type& prior)
        : variance_updater_{prior.variance.prior}
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

        std::normal_distribution<double> normal_dist;
        const auto normal_prior_for_marker
            = detail::make_normal_prior_provider<Kind>(variance);

        double sum_squares = 0.0;
        for (const Eigen::Index marker : valid_indices)
        {
            const double old_value = coefficients(marker);
            const auto likelihood = detail::make_coefficient_likelihood(
                projection, marker, old_value, residual);
            const auto posterior
                = (likelihood + normal_prior_for_marker(marker))
                      .normal_parameters();

            const double new_value = normal_dist(rng, posterior);
            if (old_value != new_value)
            {
                projection.axpy(
                    marker, old_value - new_value, residual.adjusted_response);
            }
            state.transition(marker, new_value);

            if constexpr (Kind == VarianceLayout::Pooled)
            {
                sum_squares += new_value * new_value;
            }
            else
            {
                variance_updater_.update(
                    variance(marker), 1, new_value * new_value, rng);
            }
        }
        if constexpr (Kind == VarianceLayout::Pooled)
        {
            variance_updater_.update(
                variance, valid_indices.size(), sum_squares, rng);
        }
    }

   private:
    detail::NormalVarianceConjugateUpdater variance_updater_;
};

template <VarianceLayout Kind>
[[nodiscard]] auto make_kernel(const GaussianPrior<Kind>& prior)
{
    return GaussianKernel<Kind>{prior};
}

}  // namespace gelex

#endif  // GELEX_BAYES_GENETIC_GAUSSIAN_KERNEL_H_
