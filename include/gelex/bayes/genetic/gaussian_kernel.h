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

#ifndef GELEX_BAYES_GENETIC_GAUSSIAN_KERNEL_H_
#define GELEX_BAYES_GENETIC_GAUSSIAN_KERNEL_H_

#include <Eigen/Core>
#include <random>

#include "gelex/bayes/design.h"
#include "gelex/bayes/genetic/prior.h"
#include "gelex/bayes/genetic/state.h"
#include "gelex/bayes/state.h"
#include "gelex/infra/stats/normal_sampler.h"
#include "gelex/infra/stats/scaled_inv_chi2_sampler.h"
#include "gelex/types/genetic_mode.h"

namespace gelex::detail
{

template <VarianceLayout Kind>
class GaussianKernel
{
    using Prior = GaussianPrior<Kind>;
    using State = GeneticModeState<
        GaussianState<Kind>,
        typename Prior::component_layout>;

   public:
    explicit GaussianKernel(const Prior& prior)
        : variance_sampler_{prior.variance.prior()}
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
        auto& variance = state.family_state.variance;

        // TODO(rlchen): Skip deterministic fitted-cache maintenance during
        // burn-in and rebuild it once at the sampling boundary.
        previous_adjusted_response_ = residual.adjusted_response;
        normal_.reset();
        variance_sampler_.reset();
        if constexpr (Kind == VarianceLayout::Pooled)
        {
            normal_.set_prior_var(variance);
        }

        double sum_squares = 0.0;
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
            const double new_value = normal_(
                NormalSampler<double>::Kernel{
                    .quadratic = xtx_diag(marker),
                    .linear = linear,
                    .scale = residual.variance},
                rng);
            coefficients(marker) = new_value;
            projection.axpy(
                marker, old_value - new_value, residual.adjusted_response);
            if constexpr (Kind == VarianceLayout::Pooled)
            {
                sum_squares += new_value * new_value;
            }
            else
            {
                variance(marker) = variance_sampler_(
                    {.n = 1, .sum_squares = new_value * new_value}, rng);
            }
        }

        state.component_fitted_values.col(0).noalias()
            += previous_adjusted_response_ - residual.adjusted_response;
        if constexpr (Kind == VarianceLayout::Pooled)
        {
            variance = variance_sampler_(
                {.n = static_cast<Eigen::Index>(valid_indices.size()),
                 .sum_squares = sum_squares},
                rng);
        }
    }

   private:
    ScaledInvChi2Sampler<double> variance_sampler_;
    NormalSampler<double> normal_{0.0};
    Eigen::VectorXd previous_adjusted_response_;
};

template <VarianceLayout Kind>
[[nodiscard]] auto make_mode_kernel(const GaussianPrior<Kind>& prior)
    -> GaussianKernel<Kind>
{
    return GaussianKernel<Kind>{prior};
}

}  // namespace gelex::detail

#endif  // GELEX_BAYES_GENETIC_GAUSSIAN_KERNEL_H_
