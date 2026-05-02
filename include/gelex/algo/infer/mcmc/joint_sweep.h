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

#ifndef GELEX_ALGO_INFER_MCMC_JOINT_SWEEP_H_
#define GELEX_ALGO_INFER_MCMC_JOINT_SWEEP_H_

#include <random>

#include <Eigen/Core>

#include "gelex/algo/infer/detail/genetic_binding.h"
#include "gelex/algo/infer/detail/marker_op.h"
#include "gelex/algo/infer/mcmc/kernels/concept.h"
#include "gelex/infra/stats/descriptive.h"
#include "gelex/model/bayes/effects.h"
#include "gelex/model/bayes/genotype_storage.h"
#include "gelex/model/bayes/states.h"

namespace gelex::mcmc
{

// NOLINTBEGIN(cppcoreguidelines-avoid-const-or-ref-data-members)
class GeneticJointSweep
{
   public:
    GeneticJointSweep(
        infer::detail::GeneticBlockDeps<bayes::GeneticState> first,
        infer::detail::GeneticBlockDeps<bayes::GeneticState> second,
        bayes::ResidualState& residual,
        std::mt19937_64& rng)
        : first_(first), second_(second), residual_(residual), rng_(rng)
    {
    }

    GeneticJointSweep(const GeneticJointSweep&) = delete;
    auto operator=(const GeneticJointSweep&) -> GeneticJointSweep& = delete;
    GeneticJointSweep(GeneticJointSweep&&) noexcept = default;
    auto operator=(GeneticJointSweep&&) -> GeneticJointSweep& = delete;
    ~GeneticJointSweep() = default;

    template <GeneticJointKernel Kernel>
    auto run(Kernel& kernel) -> void
    {
        const auto& first_x = bayes::get_matrix_ref(first_.effect.X);
        const auto& second_x = bayes::get_matrix_ref(second_.effect.X);
        Eigen::VectorXd& first_coeffs = first_.state.coeffs;
        Eigen::VectorXd& second_coeffs = second_.state.coeffs;
        Eigen::VectorXd& first_u = first_.state.u;
        Eigen::VectorXd& second_u = second_.state.u;
        Eigen::VectorXd& y_adj = residual_.y_adj;
        const double residual_variance = residual_.variance;

        kernel.prepare();

        for (Eigen::Index i = 0; i < first_coeffs.size(); ++i)
        {
            if (first_.effect.is_monomorphic(i)
                || second_.effect.is_monomorphic(i))
            {
                continue;
            }

            const double first_old = first_coeffs(i);
            const double second_old = second_coeffs(i);
            const auto first_col = first_x.col(i);
            const auto second_col = second_x.col(i);
            const double first_xtx = first_.effect.XtX_diag(i);
            const double second_xtx = second_.effect.XtX_diag(i);
            const double first_rhs = infer::detail::blas_ddot(first_col, y_adj)
                                     + (first_xtx * first_old);
            const double second_rhs
                = infer::detail::blas_ddot(second_col, y_adj)
                  + (second_xtx * second_old);

            const auto [first_new, second_new] = kernel.sample(
                i,
                first_xtx,
                first_rhs,
                second_xtx,
                second_rhs,
                residual_variance,
                rng_);

            first_coeffs(i) = first_new;
            second_coeffs(i) = second_new;

            infer::detail::apply_marker_update(
                y_adj,
                first_u,
                {},
                first_col,
                {.old_value = first_old, .new_value = first_new});
            infer::detail::apply_marker_update(
                y_adj,
                second_u,
                {},
                second_col,
                {.old_value = second_old, .new_value = second_new});
        }

        kernel.commit(rng_);

        first_.state.variance = gelex::stats::detail::var(first_u)(0);
        second_.state.variance = gelex::stats::detail::var(second_u)(0);
    }

   private:
    infer::detail::GeneticBlockDeps<bayes::GeneticState> first_;
    infer::detail::GeneticBlockDeps<bayes::GeneticState> second_;
    bayes::ResidualState& residual_;
    std::mt19937_64& rng_;
};
// NOLINTEND(cppcoreguidelines-avoid-const-or-ref-data-members)

}  // namespace gelex::mcmc

#endif  // GELEX_ALGO_INFER_MCMC_JOINT_SWEEP_H_
