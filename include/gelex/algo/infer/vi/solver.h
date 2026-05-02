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

#ifndef GELEX_ALGO_INFER_VI_SOLVER_H_
#define GELEX_ALGO_INFER_VI_SOLVER_H_

#include <cmath>
#include <cstddef>
#include <limits>
#include <numbers>
#include <variant>

#include <omp.h>
#include <Eigen/Core>

#include "gelex/algo/infer/params.h"
#include "gelex/algo/infer/posterior_calculator.h"
#include "gelex/infra/eigen_thread_guard.h"
#include "gelex/infra/logging/fit_event.h"
#include "gelex/infra/logging/notify.h"
#include "gelex/model/bayes/model.h"
#include "gelex/model/bayes/prior.h"
#include "gelex/types/vi_result.h"

namespace gelex
{

namespace detail
{

// Scaled-Inv-χ²(ν, s²) log density (up to additive constant)
inline auto log_scaled_inv_chi_sq(double x, double nu, double s2) -> double
{
    return -(nu / 2.0 + 1.0) * std::log(x) - (nu * s2) / (2.0 * x);
}

inline auto compute_elbo(
    const BayesModel& model,
    const bayes::Priors& priors,
    const vi::State& state) -> double
{
    const double sigma2_e = state.residual().variance;
    const auto N = static_cast<double>(model.num_individuals());

    // --- Likelihood: E_q[log p(y|β, σ²_e)] ---
    double expected_rss = state.residual().y_adj.squaredNorm();
    for (const auto& gs : state.genetics())
    {
        const auto* effect = model.genetic(gs.type);
        expected_rss += (effect->XtX_diag.array() * gs.sigma2.array()).sum();
    }
    double elbo = -0.5 * N * std::log(2.0 * std::numbers::pi * sigma2_e)
                  - 0.5 / sigma2_e * expected_rss;

    // --- Genetic prior + entropy per genetic effect ---
    for (const auto& gs : state.genetics())
    {
        const auto* prior = priors.genetic(gs.type);
        const auto& marker_prior
            = std::get<bayes::ContinuousPrior>(prior->marker);
        const double sigma2_beta = gs.marker_variance(0);
        const auto p = static_cast<double>(gs.coeffs.size());

        // E_q[log p(β|σ²_β)]
        const double sum_sq
            = (gs.coeffs.array().square() + gs.sigma2.array()).sum();
        elbo += -0.5 * p * std::log(2.0 * std::numbers::pi * sigma2_beta)
                - 0.5 / sigma2_beta * sum_sq;

        // H[q(β)] = Σ 1/2 (1 + log(2π σ²_j))
        elbo += 0.5 * gs.sigma2.array().log().sum()
                + 0.5 * p * (1.0 + std::log(2.0 * std::numbers::pi));

        // log p(σ²_β | ν, s²) — variance prior density
        elbo += log_scaled_inv_chi_sq(
            sigma2_beta,
            marker_prior.variance.param.nu,
            marker_prior.variance.param.s2);
    }

    // --- Residual variance prior ---
    elbo += log_scaled_inv_chi_sq(
        sigma2_e, priors.residual().param.nu, priors.residual().param.s2);

    return elbo;
}

}  // namespace detail

namespace vi
{

template <typename TraitUpdater>
class Solver
{
   public:
    explicit Solver(vi::Params params, TraitUpdater updater);

    auto run(
        const BayesModel& model,
        const bayes::Priors& priors,
        const FitObserver& observer = {}) -> vi::Result;

   private:
    vi::Params params_;
    TraitUpdater updater_;
};

template <typename TraitUpdater>
Solver<TraitUpdater>::Solver(vi::Params params, TraitUpdater updater)
    : params_(params), updater_(std::move(updater))
{
}

template <typename TraitUpdater>
auto Solver<TraitUpdater>::run(
    const BayesModel& model,
    const bayes::Priors& priors,
    const FitObserver& observer) -> vi::Result
{
    notify(observer, FitPriorSetEvent{.priors = &priors});

    vi::State state{model, priors};

    const ::gelex::infra::detail::EigenThreadGuard guard;
    omp_set_num_threads(1);

    double prev_elbo = -std::numeric_limits<double>::max();

    for (Eigen::Index iter = 0; iter < params_.max_iters; ++iter)
    {
        updater_(model, priors, state);
        state.compute_heritability();

        const double elbo = detail::compute_elbo(model, priors, state);
        const double delta
            = std::abs(elbo - prev_elbo) / (std::abs(prev_elbo) + 1e-300);
        const bool converged = iter > 0 && delta < params_.tol;

        notify(
            observer,
            FitVIProgressEvent{
                .current = static_cast<size_t>(iter + 1),
                .elbo = elbo,
                .delta = delta,
            });

        if (converged)
        {
            break;
        }
        prev_elbo = elbo;
    }

    notify(
        observer,
        FitVIProgressEvent{
            .done = true,
        });

    vi::Result result(std::move(state), model);
    notify(observer, FitVICompleteEvent{&result, &model});
    return result;
}

}  // namespace vi

}  // namespace gelex

#endif  // GELEX_ALGO_INFER_VI_SOLVER_H_
