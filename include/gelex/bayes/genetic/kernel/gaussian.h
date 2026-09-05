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

#ifndef GELEX_BAYES_GENETIC_KERNEL_GAUSSIAN_H_
#define GELEX_BAYES_GENETIC_KERNEL_GAUSSIAN_H_

#include <Eigen/Core>
#include <random>

#include "gelex/bayes/detail/normal_variance_conjugate_updater.h"
#include "gelex/bayes/genetic/detail/coefficient_likelihood.h"
#include "gelex/bayes/genetic/detail/normal_prior_provider.h"
#include "gelex/bayes/genetic/prior.h"
#include "gelex/bayes/genetic/state.h"
#include "gelex/bayes/genotype/design.h"
#include "gelex/bayes/state.h"
#include "gelex/genetic_mode.h"

namespace gelex::detail
{

template <VarianceLayout Kind>
class GaussianKernel
{
    using Prior = GaussianPrior<Kind>;
    using State = GeneticModeState<GaussianState<Kind>>;

   public:
    explicit GaussianKernel(const Prior& prior)
        : variance_updater_{prior.variance.prior}
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
        auto& variance = state.family_state.variance;

        // TODO(rlchen): Skip deterministic fitted-cache maintenance during
        // burn-in and rebuild it once at the sampling boundary.
        previous_adjusted_response_ = residual.adjusted_response;
        std::normal_distribution<double> normal_dist;
        const auto normal_prior_for_marker
            = make_normal_prior_provider<Kind>(variance);

        double sum_squares = 0.0;
        for (const Eigen::Index marker : valid_indices)
        {
            const double old_value = coefficients(marker);
            const auto likelihood = make_coefficient_likelihood(
                projection, marker, old_value, residual);
            const auto posterior
                = (likelihood + normal_prior_for_marker(marker))
                      .normal_parameters();

            const double new_value = normal_dist(rng, posterior);
            coefficients(marker) = new_value;
            projection.axpy(
                marker, old_value - new_value, residual.adjusted_response);

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

        state.family_state.fitted_values.col(0).noalias()
            += previous_adjusted_response_ - residual.adjusted_response;
        if constexpr (Kind == VarianceLayout::Pooled)
        {
            variance_updater_.update(
                variance, valid_indices.size(), sum_squares, rng);
        }
    }

   private:
    NormalVarianceConjugateUpdater variance_updater_;
    Eigen::VectorXd previous_adjusted_response_;
};

}  // namespace gelex::detail

#endif  // GELEX_BAYES_GENETIC_KERNEL_GAUSSIAN_H_
