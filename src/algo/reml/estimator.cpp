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

#include <Eigen/Core>
#include <cstddef>
#include <fmt/format.h>
#include <ranges>
#include <string>
#include <utility>
#include <vector>

#include "gelex/algo/reml/constrain.h"
#include "gelex/algo/reml/policy.h"
#include "gelex/algo/reml/statistics.h"
#include "gelex/algo/reml/summary.h"
#include "gelex/algo/reml/variance_calculator.h"
#include "gelex/freq/model.h"
#include "gelex/infra/logging/notify.h"
#include "gelex/infra/logging/reml_event.h"

namespace gelex
{

// Backtracking gives up below this step length, treating the anchor as
// stationary.
constexpr double LINE_SEARCH_MIN_STEP = 1e-4;

// Armijo sufficient-decrease constant c1 (Nocedal & Wright eq. 3.4).
constexpr double ARMIJO_C1 = 1e-4;

Estimator::Estimator(size_t max_iter, double tol, RemlObserver observer)
    : convergence_checker_(tol),
      max_iter_(max_iter),
      observer_(std::move(observer))
{
}

auto Estimator::fit(
    const gelex::FreqModel& model,
    gelex::FreqState& state,
    bool em_init) -> RemlFit
{
    convergence_checker_.clear();

    RemlBuffer buffer(model);

    Eigen::ArrayX<bool> constrained_mask;

    // EM initialization: one fixed-point step to seed good AI priors.
    if (em_init)
    {
        evaluate_point(model, state, buffer);
        Eigen::VectorXd sigma = EMPolicy::apply(model, state, buffer);
        constrained_mask = constrain(sigma, buffer.phenotype_variance());
        distribute_variance_components(state, sigma);
    }

    // AI-REML with Armijo backtracking: backtrack the step from 1 along the AI
    // direction p until the sufficient-increase condition holds,
    //   logL(anchor + step*p) >= logL(anchor) + c1 * step * gradᵀp,
    // so overshoot on the ill-conditioned likelihood ridge of collinear GRMs
    // cannot spiral into a limit cycle.
    Eigen::VectorXd anchor_sigma = collect_variance_components(state);
    double anchor_loglike = evaluate_point(model, state, buffer);

    double loglike = anchor_loglike;
    bool converged = false;
    size_t iter_count = max_iter_;

    for (size_t iter = 1; iter <= max_iter_; ++iter)
    {
        const Eigen::VectorXd direction = AIPolicy::direction(model, buffer);
        const double directional_derivative = buffer.first_grad.dot(direction);

        double step = 1.0;
        Eigen::VectorXd sigma;
        bool improved = false;
        while (true)
        {
            sigma = anchor_sigma + step * direction;
            constrained_mask = constrain(sigma, buffer.phenotype_variance());
            distribute_variance_components(state, sigma);
            loglike = evaluate_point(model, state, buffer);
            if (loglike
                >= anchor_loglike + (ARMIJO_C1 * step * directional_derivative))
            {
                improved = true;
                break;
            }
            if (step < LINE_SEARCH_MIN_STEP)
            {
                break;
            }
            step *= 0.5;
        }

        // No sufficient increase even at the smallest step: the anchor is
        // stationary. Restore it (refilling buffer for the final effect/SE
        // pass) and stop.
        if (!improved)
        {
            sigma = anchor_sigma;
            loglike = anchor_loglike;
            constrained_mask = constrain(sigma, buffer.phenotype_variance());
            distribute_variance_components(state, sigma);
            evaluate_point(model, state, buffer);
            converged = true;
        }

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
                .loglike = loglike,
                .labels = std::move(labels),
                .variances = std::move(variances)});

        if (!improved)
        {
            iter_count = iter;
            break;
        }

        if (convergence_checker_.is_converged(sigma, loglike))
        {
            converged = true;
            iter_count = iter;
            break;
        }

        anchor_sigma = sigma;
        anchor_loglike = loglike;
    }

    // The constraint mask (order: residual, random[0..]) marks components
    // clamped to the boundary; flag the random ones so their Wald test is
    // suppressed. The residual (index 0) is never tested against the boundary.
    for (auto&& [i, r] : std::views::enumerate(state.random()))
    {
        r.at_boundary = constrained_mask(1 + i);
    }

    const auto num_constrained = constrained_mask.count();
    const auto num_total = static_cast<Eigen::Index>(state.random().size()) + 1;
    if (2 * num_constrained > num_total)
    {
        notify(
            observer_,
            RemlConstrainedEvent{
                .num_constrained = static_cast<size_t>(num_constrained),
                .num_total = static_cast<size_t>(num_total)});
    }

    return assemble_reml_fit(
        model, state, buffer, loglike, converged, iter_count);
}

}  // namespace gelex
