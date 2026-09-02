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

#include "gelex/freq/reml/estimator.h"

#include <Eigen/Core>
#include <cstddef>
#include <fmt/format.h>
#include <string>
#include <utility>
#include <vector>

#include "gelex/freq/model.h"
#include "gelex/freq/reml/constrain.h"
#include "gelex/freq/reml/policy.h"
#include "gelex/freq/reml/progress.h"
#include "gelex/freq/reml/statistics.h"
#include "gelex/freq/reml/summary.h"
#include "gelex/freq/reml/variance_calculator.h"
#include "gelex/infra/notify.h"

namespace gelex
{

// Backtracking gives up below this step length, leaving the iteration stalled
// at the anchor.
constexpr double line_search_min_step = 1e-4;

// Armijo sufficient-decrease constant c1 (Nocedal & Wright eq. 3.4).
constexpr double armijo_c1 = 1e-4;

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

    if (em_init)
    {
        evaluate_point(model, state, buffer);
        Eigen::VectorXd sigma = EMPolicy::apply(model, state, buffer);
        constrain(sigma, buffer.phenotype_variance());
        distribute_variance_components(state, sigma);
    }

    Eigen::VectorXd anchor_sigma = collect_variance_components(state);
    double anchor_loglike = evaluate_point(model, state, buffer);

    double loglike = anchor_loglike;
    bool converged = false;
    size_t iter_count = max_iter_;

    for (size_t iter = 1; iter <= max_iter_; ++iter)
    {
        const Eigen::VectorXd direction = AIPolicy::direction(model, buffer);
        const Eigen::VectorXd anchor_grad = buffer.first_grad;

        double step = 1.0;
        Eigen::VectorXd sigma;
        bool improved = false;

        if (anchor_grad.dot(direction) > 0.0)
        {
            while (true)
            {
                sigma = anchor_sigma + step * direction;
                constrain(sigma, buffer.phenotype_variance());
                distribute_variance_components(state, sigma);
                loglike = evaluate_point(model, state, buffer);
                if (loglike >= anchor_loglike
                                   + (armijo_c1
                                      * anchor_grad.dot(sigma - anchor_sigma)))
                {
                    improved = true;
                    break;
                }
                if (step < line_search_min_step)
                {
                    break;
                }
                step *= 0.5;
            }
        }

        if (!improved)
        {
            sigma = anchor_sigma;
            loglike = anchor_loglike;
            constrain(sigma, buffer.phenotype_variance());
            distribute_variance_components(state, sigma);
            evaluate_point(model, state, buffer);
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

    const double floor_value = constraint_floor(buffer.phenotype_variance());

    auto num_constrained
        = static_cast<size_t>(state.residual().variance <= floor_value);
    for (auto& r : state.random())
    {
        r.at_boundary = r.variance <= floor_value;
        num_constrained += static_cast<size_t>(r.at_boundary);
    }

    const size_t num_total = state.random().size() + 1;
    if (2 * num_constrained > num_total)
    {
        notify(
            observer_,
            RemlConstrainedEvent{
                .num_constrained = num_constrained, .num_total = num_total});
    }

    return assemble_reml_fit(
        model, state, buffer, loglike, converged, iter_count);
}

}  // namespace gelex
