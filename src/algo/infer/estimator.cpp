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

#include "gelex/algo/infer/estimator.h"

#include <fmt/format.h>
#include <Eigen/Core>

#include "gelex/algo/infer/effect_solver.h"
#include "gelex/algo/numerics/policy.h"
#include "gelex/algo/numerics/variance_calculator.h"
#include "gelex/algo/stats/statistics.h"
#include "gelex/infra/logging/notify.h"
#include "gelex/model/freq/model.h"

namespace gelex
{

Estimator::Estimator(size_t max_iter, double tol, RemlObserver observer)
    : optimizer_(tol), max_iter_(max_iter), observer_(std::move(observer))
{
}

auto Estimator::fit(const FreqModel& model, FreqState& state, bool em_init)
    -> Eigen::MatrixXd
{
    OptimizerState opt_state(model);

    // EM initialization
    if (em_init)
    {
        em_step(model, state, opt_state);
    }

    // AI iterations
    for (size_t iter = 1; iter <= max_iter_; ++iter)
    {
        optimizer_.step<AIPolicy>(model, state, opt_state);

        loglike_ = variance_calculator::compute_loglike(model, opt_state);

        std::vector<std::string> labels;
        std::vector<double> variances;
        for (const auto& g : state.genetic())
        {
            labels.push_back(fmt::format("V({})", g.type));
            variances.push_back(g.variance);
        }
        for (const auto& r : state.random())
        {
            labels.push_back(fmt::format("V({})", r.name));
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

    // compute final results
    effect_solver::compute_fixed_effects(model, state, opt_state);
    effect_solver::compute_random_effects(model, state, opt_state);
    statistics::compute_variance_se(state, opt_state);
    statistics::compute_heritability(state, opt_state);

    notify(
        observer_,
        RemlCompleteEvent{
            .model = &model,
            .state = &state,
            .converged = converged_,
            .iter_count = iter_count_,
            .max_iter = max_iter_,
            .loglike = loglike_});

    return std::move(opt_state.v);
}

auto Estimator::em_step(
    const FreqModel& model,
    FreqState& state,
    OptimizerState& opt_state) -> void
{
    optimizer_.step<EMPolicy>(model, state, opt_state);

    double loglike = variance_calculator::compute_loglike(model, opt_state);

    std::vector<double> init_variances;
    for (const auto& g : state.genetic())
    {
        init_variances.push_back(g.variance);
    }
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
