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

#include <array>
#include <cstdint>
#include <random>
#include <vector>

#include <Eigen/Core>

#include "gelex/algo/infer/detail/marker_op.h"
#include "gelex/algo/infer/mcmc/kernels/common.h"
#include "gelex/algo/infer/mcmc/kernels/mixture_op.h"
#include "gelex/infra/stats/conjugate_prior.h"
#include "gelex/model/bayes/effects.h"
#include "gelex/model/bayes/genotype_storage.h"
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
        const bayes::GeneticEffect& effect)
        : state_(state),
          assignment_(
              unpack_marker_allocation<bayes::ComponentAllocation>(
                  state,
                  "BayesRKernel")),
          X_(bayes::get_matrix_ref(effect.X)),
          multiplier_(
              unpack_marker_prior<bayes::MixturePrior>(prior, "BayesRKernel")
                  .multiplier),
          marker_variances_(multiplier_.size()),
          logpi_(multiplier_.size()),
          pi_count_(Eigen::VectorXi::Zero(multiplier_.size())),
          variance_sampler_(
              unpack_marker_prior<bayes::MixturePrior>(prior, "BayesRKernel")
                  .variance.param.nu,
              unpack_marker_prior<bayes::MixturePrior>(prior, "BayesRKernel")
                  .variance.param.s2),
          monomorphic_indices_(collect_monomorphic_indices(effect))
    {
    }

    auto prepare() -> void
    {
        variance_sampler_.reset();
        uniform_.reset();
        normal_.reset();

        logpi_ = assignment_.assignment.proportion.array().log();
        marker_variances_ = state_.marker_variance(0) * multiplier_.array();

        pi_count_.setZero();
        for (Eigen::Index i : monomorphic_indices_)
        {
            pi_count_(assignment_.assignment.tracker(i))++;
        }
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
        const double old_i = state_.coeffs(marker_index);

        const auto scale_pp
            = compute_scale_posteriors(rhs, xtx_diag, residual_variance);
        const int8_t dist_index
            = draw_component(scale_pp, num_scale_classes, rng);

        // NOLINTBEGIN(cppcoreguidelines-pro-bounds-constant-array-index)
        const double new_i = (dist_index > 0)
                                 ? sample_slab_effect(scale_pp[dist_index], rng)
                                 : 0.0;
        // NOLINTEND(cppcoreguidelines-pro-bounds-constant-array-index)

        update_accumulator(dist_index, new_i, marker_index);

        // NOLINTNEXTLINE(readability-suspicious-call-argument)
        detail::update_component_u(
            assignment_.component_u,
            old_index,
            old_i,
            dist_index,
            new_i,
            X_.col(marker_index));

        return new_i;
    }

    auto commit(std::mt19937_64& rng) -> void
    {
        const Eigen::Index p = state_.coeffs.size();
        assignment_.assignment.count = pi_count_;
        const Eigen::Index num_nonzero = p - assignment_.assignment.count(0);

        state_.marker_variance(0)
            = variance_sampler_({num_nonzero, sum_square_coeffs_}, rng);

        detail::compute_component_variances(assignment_);
    }

   private:
    bayes::GeneticState& state_;
    bayes::ComponentAllocation& assignment_;
    Eigen::Ref<const Eigen::MatrixXd> X_;
    Eigen::VectorXd multiplier_;
    Eigen::VectorXd marker_variances_;
    Eigen::VectorXd logpi_;
    Eigen::VectorXi pi_count_;
    double sum_square_coeffs_{0.0};

    ScaledInvChi2Sampler<double> variance_sampler_;
    std::uniform_real_distribution<double> uniform_{0.0, 1.0};
    std::normal_distribution<double> normal_{0.0, 1.0};

    std::vector<Eigen::Index> monomorphic_indices_;

    static auto collect_monomorphic_indices(const bayes::GeneticEffect& effect)
        -> std::vector<Eigen::Index>
    {
        const Eigen::Index p = bayes::get_cols(effect.X);
        std::vector<Eigen::Index> indices;
        for (Eigen::Index i = 0; i < p; ++i)
        {
            if (effect.is_monomorphic(i))
            {
                indices.push_back(i);
            }
        }
        return indices;
    }

    auto compute_scale_posteriors(
        double rhs,
        double xtx_diag,
        double residual_variance) const
        -> std::array<detail::PosteriorParams, detail::kMaxMixtureComponents>
    {
        std::array<detail::PosteriorParams, detail::kMaxMixtureComponents>
            scale_pp{};
        const Eigen::Index num_scale_classes = multiplier_.size();
        // NOLINTBEGIN(cppcoreguidelines-pro-bounds-constant-array-index,readability-identifier-length)
        for (Eigen::Index c = 1; c < num_scale_classes; ++c)
        {
            scale_pp[c] = detail::compute_posterior_params(
                rhs, marker_variances_(c), xtx_diag, residual_variance);
        }
        // NOLINTEND(cppcoreguidelines-pro-bounds-constant-array-index,readability-identifier-length)
        return scale_pp;
    }

    auto draw_component(
        const std::array<
            detail::PosteriorParams,
            detail::kMaxMixtureComponents>& scale_pp,
        Eigen::Index num_scale_classes,
        std::mt19937_64& rng) -> int8_t
    {
        // NOLINTBEGIN(cppcoreguidelines-pro-bounds-constant-array-index,readability-identifier-length)
        std::array<double, detail::kMaxMixtureComponents> ll{};
        ll[0] = logpi_(0);
        double max_ll = ll[0];
        for (Eigen::Index c = 1; c < num_scale_classes; ++c)
        {
            ll[c] = scale_pp[c].log_likelihood_kernel + logpi_(c);
            max_ll = std::max(ll[c], max_ll);
        }

        std::array<double, detail::kMaxMixtureComponents> probs{};
        double total = 0.0;
        for (Eigen::Index k = 0; k < num_scale_classes; ++k)
        {
            probs[k] = std::exp(ll[k] - max_ll);
            total += probs[k];
        }

        const double u_val = uniform_(rng) * total;
        auto dist_index = static_cast<int8_t>(num_scale_classes - 1);
        double cumsum = 0.0;
        for (Eigen::Index k = 0; k < num_scale_classes; ++k)
        {
            cumsum += probs[k];
            if (u_val < cumsum)
            {
                dist_index = static_cast<int8_t>(k);
                break;
            }
        }
        // NOLINTEND(cppcoreguidelines-pro-bounds-constant-array-index,readability-identifier-length)
        return dist_index;
    }

    auto sample_slab_effect(
        const detail::PosteriorParams& pp,
        std::mt19937_64& rng) -> double
    {
        return (normal_(rng) * pp.stddev) + pp.mean;
    }

    auto update_accumulator(
        int8_t dist_index,
        double new_i,
        Eigen::Index marker_index) -> void
    {
        assignment_.assignment.tracker(marker_index) = dist_index;
        pi_count_(dist_index)++;
        if (dist_index > 0)
        {
            sum_square_coeffs_ += (new_i * new_i) / multiplier_(dist_index);
        }
    }
};
// NOLINTEND(cppcoreguidelines-avoid-const-or-ref-data-members)

}  // namespace gelex::mcmc

#endif  // GELEX_ALGO_INFER_MCMC_KERNELS_BAYES_R_H_
