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

#ifndef GELEX_ALGO_INFER_MCMC_KERNELS_BAYES_A_H_
#define GELEX_ALGO_INFER_MCMC_KERNELS_BAYES_A_H_

#include <random>
#include <variant>

#include "gelex/algo/infer/mcmc/kernels/common.h"
#include "gelex/algo/infer/mcmc/state.h"
#include "gelex/infra/stats/conjugate_prior.h"

namespace gelex::mcmc
{

// BayesA: per-marker scaled-inv-chi^2 variance, no spike.
// Mutates: state.coeffs (via sweep), state.marker_variance(i).
// NOLINTBEGIN(cppcoreguidelines-avoid-const-or-ref-data-members)
class BayesAKernel
{
   public:
    BayesAKernel(const bayes::GeneticPrior& prior, bayes::GeneticState& state)
        : state_(state),
          normal_(0.0),
          sigma_(make_variance_sampler(
              std::get<bayes::GeneticSpec>(prior.spec).variance))
    {
    }

    auto prepare() -> void
    {
        normal_.reset();
        sigma_.reset();
    }

    auto sample(
        Eigen::Index marker_index,
        double xtx_diag,
        double rhs,
        double residual_variance,
        std::mt19937_64& rng) -> double
    {
        normal_.set_prior_var(state_.marker_variance(marker_index));
        const double new_i = normal_(
            stats::NormalSampler<double>::Kernel{
                .quadratic = xtx_diag,
                .linear = rhs,
                .scale = residual_variance,
            },
            rng);
        state_.marker_variance(marker_index) = sigma_({1, new_i * new_i}, rng);
        return new_i;
    }

    auto commit(std::mt19937_64& /*rng*/) -> void {}

   private:
    bayes::GeneticState& state_;
    stats::NormalSampler<double> normal_;
    stats::ScaledInvChi2Sampler<double> sigma_;
};
// NOLINTEND(cppcoreguidelines-avoid-const-or-ref-data-members)

}  // namespace gelex::mcmc

#endif  // GELEX_ALGO_INFER_MCMC_KERNELS_BAYES_A_H_
