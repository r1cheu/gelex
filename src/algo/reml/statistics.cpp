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

#include "gelex/algo/reml/statistics.h"

#include <algorithm>
#include <cmath>
#include <cstddef>

#include <Eigen/Core>

#include "gelex/algo/reml/optimizer_state.h"
#include "gelex/freq/model.h"

namespace gelex::reml
{

auto compute_aic(const FreqModel& model, double loglike) -> double
{
    // k = number of variance components + number of fixed effects
    auto k = static_cast<double>(
        1 + model.random().size() + model.fixed().X.cols());
    return -2.0 * loglike + 2.0 * k;
}

auto compute_bic(const FreqModel& model, double loglike) -> double
{
    auto k = static_cast<double>(
        1 + model.random().size() + model.fixed().X.cols());
    auto n = static_cast<double>(model.num_individuals());
    return -2.0 * loglike + k * std::log(n);
}

auto compute_variance_se(FreqState& state, const OptimizerState& opt_state)
    -> void
{
    // se(σ) = sqrt(diag(-H⁻¹))
    // variance component order: residual, random[0..]
    Eigen::VectorXd se = (-opt_state.hess_inv.diagonal()).array().sqrt();

    Eigen::Index idx = 0;

    // residual
    state.residual().variance_se = se(idx++);

    // random effects
    for (auto& r : state.random())
    {
        r.variance_se = se(idx++);
    }
}

auto compute_variance_ratio(FreqState& state, const OptimizerState& opt_state)
    -> void
{
    // total phenotypic variance
    double sum_var = state.residual().variance;
    for (const auto& r : state.random())
    {
        sum_var += r.variance;
    }

    if (sum_var <= 0.0)
    {
        return;
    }

    double sum_var_sq = sum_var * sum_var;
    auto n_comp = opt_state.hess_inv.rows();

    for (size_t ri = 0; ri < state.random().size(); ++ri)
    {
        auto& r = state.random()[ri];
        Eigen::Index r_idx = 1 + static_cast<Eigen::Index>(ri);

        r.variance_ratio = r.variance / sum_var;

        // gradient for delta method:
        // ∂ratio/∂σ_i = -σ_r / (Σσ)²           for i ≠ r
        // ∂ratio/∂σ_r = (Σσ - σ_r) / (Σσ)²     for i == r
        Eigen::VectorXd grad(n_comp);
        for (Eigen::Index i = 0; i < n_comp; ++i)
        {
            if (i == r_idx)
            {
                grad(i) = (sum_var - r.variance) / sum_var_sq;
            }
            else
            {
                grad(i) = -r.variance / sum_var_sq;
            }
        }

        double var_ratio = grad.dot(-opt_state.hess_inv * grad);
        r.variance_ratio_se = std::sqrt(std::max(0.0, var_ratio));
    }
}

}  // namespace gelex::reml
