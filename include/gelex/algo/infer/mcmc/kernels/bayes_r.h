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

#ifndef GELEX_ALGO_INFER_MCMC_KERNELS_BAYES_R_H_
#define GELEX_ALGO_INFER_MCMC_KERNELS_BAYES_R_H_

#include <cstdint>
#include <random>

#include <Eigen/Core>

#include "gelex/algo/infer/mcmc/kernels/common.h"
#include "gelex/algo/infer/mcmc/kernels/mixture_op.h"
#include "gelex/infra/stats/conjugate_prior.h"
#include "gelex/model/bayes/effects.h"
#include "gelex/model/bayes/prior.h"
#include "gelex/model/bayes/states.h"

namespace gelex::mcmc
{

// BayesR: finite mixture of scaled normals (symmetric sampler, no AT policy).
// Mutates:
//   state.coeffs (via sweep)
//   state.marker_variance(0)
//   state.group->ComponentAllocation::{tracker, count, component_u,
//                                      component_variance}
//   state.variance (via sweep)
// NOLINTBEGIN(cppcoreguidelines-avoid-const-or-ref-data-members)
class BayesRKernel
{
   public:
    BayesRKernel(
        const bayes::GeneticPrior& prior,
        bayes::GeneticState& state,
        const bayes::GeneticEffect& /*effect*/)
        : state_(state),
          assignment_(
              unpack_marker_allocation<bayes::ComponentAllocation>(
                  state,
                  "BayesRKernel")),
          multiplier_(
              unpack_marker_prior<bayes::MixturePrior>(prior, "BayesRKernel")
                  .multiplier),
          marker_variances_(multiplier_.size()),
          logpi_(multiplier_.size()),
          pi_count_(assignment_.assignment.count),
          variance_sampler_(
              make_variance_sampler<bayes::MixturePrior>(
                  prior,
                  "BayesRKernel")),
          normal_(0.0)
    {
    }

    auto prepare() -> void
    {
        variance_sampler_.reset();
        uniform_.reset();
        normal_.reset();

        logpi_ = assignment_.assignment.proportion.array().log();
        marker_variances_ = state_.marker_variance(0) * multiplier_.array();
        sum_square_coeffs_ = 0.0;
    }

    auto sample(
        Eigen::Index marker_index,
        double xtx_diag,
        double rhs,
        double residual_variance,
        std::mt19937_64& rng) -> double
    {
        const Eigen::Index num_scale_classes = multiplier_.size();
        const int8_t old_index = assignment_.assignment.tracker(marker_index);

        compute_scale_posteriors(rhs, xtx_diag, residual_variance);
        const int8_t dist_index = draw_component(num_scale_classes, rng);

        const double new_i = (dist_index > 0)
                                 ? normal_.draw(
                                       {.mean = scale_pp_.means(dist_index),
                                        .var = scale_pp_.vars(dist_index)},
                                       rng)
                                 : 0.0;

        assignment_.assignment.tracker(marker_index) = dist_index;
        pi_count_(old_index)--;
        pi_count_(dist_index)++;
        if (dist_index > 0)
        {
            sum_square_coeffs_ += (new_i * new_i) / multiplier_(dist_index);
        }
        return new_i;
    }

    auto commit(std::mt19937_64& rng) -> void
    {
        const Eigen::Index p = state_.coeffs.size();
        assignment_.assignment.count = pi_count_;
        const Eigen::Index num_nonzero = p - assignment_.assignment.count(0);

        state_.marker_variance(0)
            = variance_sampler_({num_nonzero, sum_square_coeffs_}, rng);
    }

   private:
    bayes::GeneticState& state_;
    bayes::ComponentAllocation& assignment_;
    Eigen::VectorXd multiplier_;
    Eigen::VectorXd marker_variances_;
    Eigen::VectorXd logpi_;
    Eigen::VectorXi pi_count_;
    double sum_square_coeffs_{0.0};

    ScaledInvChi2Sampler<double> variance_sampler_;
    std::uniform_real_distribution<double> uniform_{0.0, 1.0};
    NormalSampler<double> normal_;
    detail::MixtureNormalPosteriors scale_pp_;

    auto compute_scale_posteriors(
        double rhs,
        double xtx_diag,
        double residual_variance) -> void
    {
        const Eigen::Index num_scale_classes = multiplier_.size();
        scale_pp_.log_likelihoods(0) = 0.0;
        for (Eigen::Index cls = 1; cls < num_scale_classes; ++cls)
        {
            const auto post = normal_.set_prior_var(marker_variances_(cls))
                                  .posterior_with_logL(
                                      {.quadratic = xtx_diag,
                                       .linear = rhs,
                                       .scale = residual_variance});
            scale_pp_.means(cls) = post.params.mean;
            scale_pp_.vars(cls) = post.params.var;
            scale_pp_.log_likelihoods(cls) = post.log_likelihood_kernel;
        }
    }

    auto draw_component(Eigen::Index num_scale_classes, std::mt19937_64& rng)
        -> int8_t
    {
        Eigen::Array<double, detail::kMaxMixtureComponents, 1> ll;
        Eigen::Array<double, detail::kMaxMixtureComponents, 1> probs;

        ll.head(num_scale_classes)
            = scale_pp_.log_likelihoods.head(num_scale_classes)
              + logpi_.head(num_scale_classes).array();
        const double max_ll = ll.head(num_scale_classes).maxCoeff();
        probs.head(num_scale_classes)
            = (ll.head(num_scale_classes) - max_ll).exp();
        const double total = probs.head(num_scale_classes).sum();

        const double threshold = uniform_(rng) * total;
        auto dist_index = static_cast<int8_t>(num_scale_classes - 1);
        double cumsum = 0.0;
        for (Eigen::Index k = 0; k < num_scale_classes; ++k)
        {
            cumsum += probs(k);
            if (threshold < cumsum)
            {
                dist_index = static_cast<int8_t>(k);
                break;
            }
        }
        return dist_index;
    }
};
// NOLINTEND(cppcoreguidelines-avoid-const-or-ref-data-members)

}  // namespace gelex::mcmc

#endif  // GELEX_ALGO_INFER_MCMC_KERNELS_BAYES_R_H_
