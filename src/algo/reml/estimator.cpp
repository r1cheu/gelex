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

#include "gelex/algo/reml/estimator.h"

#include <fmt/format.h>
#include <Eigen/Core>
#include <cstddef>
#include <string>
#include <utility>
#include <vector>

#include "gelex/algo/reml/effect_solver.h"
#include "gelex/algo/reml/policy.h"
#include "gelex/algo/reml/result.h"
#include "gelex/algo/reml/statistics.h"
#include "gelex/algo/reml/variance_calculator.h"
#include "gelex/freq/model.h"
#include "gelex/infra/logging/notify.h"
#include "gelex/infra/logging/reml_event.h"

namespace gelex
{

Estimator::Estimator(size_t max_iter, double tol, RemlObserver observer)
    : optimizer_(tol), max_iter_(max_iter), observer_(std::move(observer))
{
}

auto Estimator::fit(
    const gelex::FreqModel& model,
    gelex::FreqState& state,
    bool em_init) -> RemlResult
{
    optimizer_.reset();
    iter_count_ = 0;
    loglike_ = 0.0;
    converged_ = false;

    OptimizerState opt_state(model);

    // EM initialization
    if (em_init)
    {
        em_step(model, state, opt_state);
        optimizer_.reset();
    }

    // AI iterations
    for (size_t iter = 1; iter <= max_iter_; ++iter)
    {
        optimizer_.step<AIPolicy>(model, state, opt_state);

        loglike_ = compute_loglike(model, opt_state);

        std::vector<std::string> labels;
        std::vector<double> variances;
        for (size_t i = 0; i < state.random().size(); ++i)
        {
            const auto& r = state.random()[i];
            labels.push_back(fmt::format("V({})", model.random()[i].name));
            variances.push_back(r.variance);
        }
        labels.emplace_back("V(e)");
        variances.push_back(state.residual().variance);

        notify(
            observer_,
            RemlIterationEvent{
                .iter = iter,
                .loglike = loglike_,
                .labels = std::move(labels),
                .variances = std::move(variances)});

        if (optimizer_.is_converged())
        {
            converged_ = true;
            iter_count_ = iter;
            break;
        }
    }

    if (!converged_)
    {
        iter_count_ = max_iter_;
    }

    const auto num_total = static_cast<Eigen::Index>(state.random().size()) + 1;
    if (2 * optimizer_.num_constrained() > num_total)
    {
        notify(
            observer_,
            RemlConstrainedEvent{
                .num_constrained
                = static_cast<size_t>(optimizer_.num_constrained()),
                .num_total = static_cast<size_t>(num_total)});
    }

    // compute final results
    compute_fixed_effects(model, state, opt_state);
    compute_random_effects(model, state, opt_state);
    compute_variance_se(state, opt_state);
    compute_variance_ratio(state, opt_state);

    // Materialize P = V^{-1} - ViX * XtViX_inv * ViX' in opt_state.V's memory
    // so downstream GWAS can use a single dense GEMM per SNP chunk.
    opt_state.V.noalias()
        -= opt_state.ViX * opt_state.XtViX_inv * opt_state.ViX.transpose();

    return RemlResult{
        .P = std::move(opt_state.V),
        .Py = std::move(opt_state.Py),
        .Vp = state.Vp()};
}

auto Estimator::em_step(
    const gelex::FreqModel& model,
    gelex::FreqState& state,
    OptimizerState& opt_state) -> void
{
    optimizer_.step<EMPolicy>(model, state, opt_state);

    double loglike = compute_loglike(model, opt_state);

    std::vector<double> init_variances;
    for (const auto& r : state.random())
    {
        init_variances.push_back(r.variance);
    }
    init_variances.push_back(state.residual().variance);

    notify(
        observer_,
        RemlEmInitEvent{
            .loglike = loglike, .init_variances = std::move(init_variances)});
}

}  // namespace gelex
