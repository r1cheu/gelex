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

#ifndef GELEX_ALGO_INFER_MCMC_KERNELS_BAYES_B_H_
#define GELEX_ALGO_INFER_MCMC_KERNELS_BAYES_B_H_

#include <cmath>
#include <cstdint>
#include <random>

#include "gelex/algo/infer/mcmc/kernels/common.h"
#include "gelex/infra/stats/conjugate_prior.h"
#include "gelex/model/bayes/prior.h"
#include "gelex/model/bayes/states.h"

namespace gelex::mcmc
{

// BayesB: spike-and-slab with per-marker scaled-inv-chi^2 slab variance.
// Mutates: state.coeffs (via sweep), state.marker_variance(i),
//          state.group->Assignment::tracker / count.
// NOLINTBEGIN(cppcoreguidelines-avoid-const-or-ref-data-members)
class BayesBKernel
{
   public:
    BayesBKernel(const bayes::GeneticPrior& prior, bayes::GeneticState& state)
        : state_(state),
          assignment_(
              unpack_marker_allocation<bayes::Assignment>(
                  state,
                  "BayesBKernel")),
          normal_(0.0),
          variance_sampler_(
              make_variance_sampler<bayes::SpikePrior>(prior, "BayesBKernel"))
    {
    }

    auto prepare() -> void
    {
        normal_.reset();
        variance_sampler_.reset();
        uniform_.reset();
        count_1_ = assignment_.count(1);
        logpi_ = assignment_.proportion.array().log();
    }

    auto sample(
        Eigen::Index marker_index,
        double xtx_diag,
        double rhs,
        double residual_variance,
        std::mt19937_64& rng) -> double
    {
        const std::int8_t old_component = assignment_.tracker(marker_index);

        const auto post
            = normal_.set_prior_var(state_.marker_variance(marker_index))
                  .posterior_with_logL(
                      stats::NormalSampler<double>::Kernel{
                          .quadratic = xtx_diag,
                          .linear = rhs,
                          .scale = residual_variance,
                      });

        const double log_like_1_minus_0
            = post.log_likelihood_kernel + logpi_(1) - logpi_(0);
        const double prob_component_0
            = 1.0 / (1.0 + std::exp(log_like_1_minus_0));

        const std::int8_t dist_index
            = (uniform_(rng) < prob_component_0) ? 0 : 1;
        assignment_.tracker(marker_index) = dist_index;
        count_1_
            += static_cast<int>(dist_index) - static_cast<int>(old_component);

        double new_i = 0.0;
        if (dist_index == 1)
        {
            new_i = normal_.draw(post.params, rng);
            state_.marker_variance(marker_index)
                = variance_sampler_({1, new_i * new_i}, rng);
        }
        return new_i;
    }

    auto commit(std::mt19937_64& /*rng*/) -> void
    {
        const auto p = static_cast<int>(state_.coeffs.size());
        assignment_.count(1) = count_1_;
        assignment_.count(0) = p - count_1_;
    }

   private:
    bayes::GeneticState& state_;
    bayes::Assignment& assignment_;
    stats::NormalSampler<double> normal_;
    stats::ScaledInvChi2Sampler<double> variance_sampler_;
    std::uniform_real_distribution<double> uniform_{0.0, 1.0};

    int count_1_{0};
    Eigen::VectorXd logpi_;
};
// NOLINTEND(cppcoreguidelines-avoid-const-or-ref-data-members)

}  // namespace gelex::mcmc

#endif  // GELEX_ALGO_INFER_MCMC_KERNELS_BAYES_B_H_
