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

#ifndef GELEX_ALGO_INFER_MCMC_KERNELS_BAYES_RR_H_
#define GELEX_ALGO_INFER_MCMC_KERNELS_BAYES_RR_H_

#include <random>

#include "gelex/algo/infer/mcmc/kernels/common.h"
#include "gelex/infra/stats/conjugate_prior.h"
#include "gelex/model/bayes/prior.h"
#include "gelex/model/bayes/states.h"

namespace gelex::mcmc
{

// BayesRR: ridge — single shared marker variance, no spike, no per-marker var.
// Mutates: state.coeffs (via sweep), state.marker_variance(0).
// NOLINTBEGIN(cppcoreguidelines-avoid-const-or-ref-data-members)
class BayesRRKernel
{
   public:
    BayesRRKernel(const bayes::GeneticPrior& prior, bayes::GeneticState& state)
        : state_(state),
          normal_(state.marker_variance(0)),
          variance_sampler_(
              make_variance_sampler<bayes::ContinuousPrior>(
                  prior,
                  "BayesRRKernel"))
    {
    }

    auto prepare() -> void
    {
        normal_.reset();
        variance_sampler_.reset();
        n_used_ = 0;
        normal_.set_prior_var(state_.marker_variance(0));
    }

    auto sample(
        Eigen::Index /*marker_index*/,
        double xtx_diag,
        double rhs,
        double residual_variance,
        std::mt19937_64& rng) -> double
    {
        const double new_i = normal_(
            NormalSampler<double>::Kernel{
                .quadratic = xtx_diag,
                .linear = rhs,
                .scale = residual_variance,
            },
            rng);
        ++n_used_;
        return new_i;
    }

    auto commit(std::mt19937_64& rng) -> void
    {
        state_.marker_variance(0)
            = variance_sampler_({n_used_, state_.coeffs.squaredNorm()}, rng);
    }

   private:
    bayes::GeneticState& state_;
    NormalSampler<double> normal_;
    ScaledInvChi2Sampler<double> variance_sampler_;
    Eigen::Index n_used_{0};
};
// NOLINTEND(cppcoreguidelines-avoid-const-or-ref-data-members)

}  // namespace gelex::mcmc

#endif  // GELEX_ALGO_INFER_MCMC_KERNELS_BAYES_RR_H_
