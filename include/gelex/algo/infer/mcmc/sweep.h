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

#include <random>

#include <Eigen/Core>

#include "gelex/algo/infer/detail/marker_op.h"
#include "gelex/algo/infer/mcmc/kernels/concept.h"
#include "gelex/infra/stats/descriptive.h"
#include "gelex/model/bayes/effects.h"
#include "gelex/model/bayes/genotype_storage.h"
#include "gelex/model/bayes/states.h"

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
        const auto& X = bayes::get_matrix_ref(effect_.X);
        const auto& XtX_diag = effect_.XtX_diag;
        Eigen::VectorXd& coeffs = state_.coeffs;
        Eigen::VectorXd& u = state_.u;
        Eigen::VectorXd& y_adj = residual_.y_adj;
        const double residual_variance = residual_.variance;

        kernel.prepare();

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
                = ::gelex::detail::blas_ddot(col, y_adj) + (xtx_diag_i * old_i);

            const double new_i
                = kernel.sample(i, xtx_diag_i, rhs, residual_variance, rng_);

            coeffs(i) = new_i;
            if (old_i != new_i)
            {
                ::gelex::detail::update_residual_and_gebv(
                    y_adj, u, col, old_i, new_i);
            }
        }

        kernel.commit(rng_);

        state_.variance = gelex::detail::var(state_.u)(0);
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
