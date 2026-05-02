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

#include "gelex/algo/infer/vi/updaters/common.h"

#include <Eigen/Core>

#include "gelex/infra/stats/conjugate_prior.h"
#include "gelex/model/bayes/effects.h"
#include "gelex/model/bayes/model.h"
#include "gelex/model/bayes/states.h"

namespace gelex::vi::detail
{

using Eigen::Index;
using Eigen::VectorXd;

auto Fixed::operator()(
    const BayesModel& model,
    const bayes::Priors& /*priors*/,
    vi::State& states) const -> void
{
    const auto& effect = model.fixed();
    auto& state = states.fixed();
    auto& residual = states.residual();

    auto& y_adj = residual.y_adj;
    auto& coeffs = state.coeffs;
    const auto& XtX_diag = effect.XtX_diag;
    const auto& X = effect.X;

    for (Index i = 0; i < coeffs.size(); ++i)
    {
        const double old_i = coeffs(i);
        const auto& col = X.col(i);
        const double norm = XtX_diag(i);

        const double rhs = col.dot(y_adj) + (norm * old_i);
        const double post_mean = rhs / norm;

        coeffs(i) = post_mean;

        const double diff = old_i - post_mean;
        y_adj.array() += diff * col.array();
    }
}

auto Random::operator()(
    const BayesModel& model,
    const bayes::Priors& priors,
    vi::State& states) const -> void
{
    if (model.random().empty())
    {
        return;
    }

    const auto& effects = model.random();
    auto& state = states.random();
    auto& residual = states.residual();

    for (Index i = 0; i < static_cast<Index>(effects.size()); ++i)
    {
        update_impl(effects[i], priors.random()[i], state[i], residual);
    }
}

auto Random::update_impl(
    const bayes::RandomEffect& effect,
    const bayes::RandomPrior& prior,
    bayes::RandomState& state,
    bayes::ResidualState& residual) -> void
{
    VectorXd& coeffs = state.coeffs;
    const VectorXd& XtX_diag = effect.XtX_diag;
    const auto& X = effect.X;

    auto& y_adj = residual.y_adj;
    const double residual_variance = residual.variance;
    const double sigma = state.variance;

    const VectorXd inv_scaler
        = 1.0 / (XtX_diag.array() + residual_variance / sigma);

    for (Index i = 0; i < coeffs.size(); ++i)
    {
        const double old_i = coeffs(i);
        const auto& col = X.col(i);
        const double norm = XtX_diag(i);

        const double rhs = col.dot(y_adj) + (norm * old_i);
        const double post_mean = rhs * inv_scaler(i);

        coeffs(i) = post_mean;

        const double diff = old_i - post_mean;
        y_adj.array() += col.array() * diff;
    }

    gelex::stats::detail::ScaledInvChiSq chi_squared{prior.param};
    chi_squared.compute(coeffs.squaredNorm(), coeffs.size());
    state.variance = chi_squared.expected_value();
}

auto Residual::operator()(
    const BayesModel& model,
    const bayes::Priors& priors,
    vi::State& states) const -> void
{
    auto& residual = states.residual();

    // E_q[||y - Xβ||²] = ||y_adj||² + tr(X'X · diag(σ²_q))
    double expected_rss = residual.y_adj.squaredNorm();
    for (const auto& gs : states.genetics())
    {
        const auto* effect = model.genetic(gs.type);
        expected_rss += (effect->XtX_diag.array() * gs.sigma2.array()).sum();
    }

    gelex::stats::detail::ScaledInvChiSq chi_squared{priors.residual().param};
    chi_squared.compute(expected_rss, model.num_individuals());
    residual.variance = chi_squared.expected_value();
}

}  // namespace gelex::vi::detail
