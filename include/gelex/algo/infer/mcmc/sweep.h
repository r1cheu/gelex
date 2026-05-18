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

#ifndef GELEX_ALGO_INFER_MCMC_SWEEP_H_
#define GELEX_ALGO_INFER_MCMC_SWEEP_H_

#include <cstdint>
#include <random>
#include <span>

#include <Eigen/Core>

#include "gelex/algo/infer/detail/marker_op.h"
#include "gelex/algo/infer/mcmc/kernels/concept.h"
#include "gelex/algo/infer/mcmc/kernels/detail/mixture_op.h"
#include "gelex/algo/infer/mcmc/state.h"
#include "gelex/data/genotype/genotype.h"
#include "gelex/infra/stats/detail/var.h"
#include "gelex/model/bayes/effects.h"

namespace gelex::mcmc
{

// NOLINTBEGIN(cppcoreguidelines-avoid-const-or-ref-data-members)
class GeneticSweep
{
   public:
    GeneticSweep(
        const bayes::GeneticEffect& effect,
        bayes::GeneticState& state,
        bayes::ResidualState& residual,
        std::mt19937_64& rng)
        : effect_(effect), state_(state), residual_(residual), rng_(rng)
    {
    }

    GeneticSweep(const GeneticSweep&) = delete;
    auto operator=(const GeneticSweep&) -> GeneticSweep& = delete;
    GeneticSweep(GeneticSweep&&) noexcept = default;
    auto operator=(GeneticSweep&&) -> GeneticSweep& = delete;
    ~GeneticSweep() = default;

    template <GeneticKernel Kernel>
    auto run(Kernel& kernel) -> void
    {
        const auto X = effect_.X.matrix();
        const auto& XtX_diag = effect_.XtX_diag;
        Eigen::VectorXd& coeffs = state_.coeffs;
        Eigen::VectorXd& u = state_.u;
        Eigen::VectorXd& y_adj = residual_.y_adj;
        const double residual_variance = residual_.variance;

        kernel.prepare();

        auto* mix = detail::get_mixture(state_);
        const bool track_components = mix && !mix->component_u.empty();
        const auto component_u_span
            = track_components ? std::span<Eigen::VectorXd>(mix->component_u)
                               : std::span<Eigen::VectorXd>{};

        for (Eigen::Index i = 0; i < coeffs.size(); ++i)
        {
            const double xtx_diag_i = XtX_diag(i);
            if (effect_.is_monomorphic(i))
            {
                continue;
            }

            const double old_i = coeffs(i);
            const auto col = X.col(i);
            const double rhs
                = infer::detail::blas_ddot(col, y_adj) + (xtx_diag_i * old_i);

            const int8_t old_class = mix ? mix->assignment.tracker(i) : -1;

            const double new_i
                = kernel.sample(i, xtx_diag_i, rhs, residual_variance, rng_);

            coeffs(i) = new_i;
            const int8_t new_class = mix ? mix->assignment.tracker(i) : -1;

            infer::detail::apply_marker_update(
                y_adj,
                u,
                component_u_span,
                col,
                {.old_value = old_i,
                 .new_value = new_i,
                 .old_class = old_class,
                 .new_class = new_class});
        }

        kernel.commit(rng_);

        state_.variance = gelex::stats::detail::var(state_.u)(0);
        if (track_components)
        {
            detail::compute_component_variances(*mix);
        }
    }

   private:
    const bayes::GeneticEffect& effect_;
    bayes::GeneticState& state_;
    bayes::ResidualState& residual_;
    std::mt19937_64& rng_;
};
// NOLINTEND(cppcoreguidelines-avoid-const-or-ref-data-members)

}  // namespace gelex::mcmc

#endif  // GELEX_ALGO_INFER_MCMC_SWEEP_H_
