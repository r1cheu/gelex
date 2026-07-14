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

#include <Eigen/Core>
#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <numbers>
#include <utility>
#include <vector>

#include "gelex/algo/reml/constrain.h"
#include "gelex/algo/reml/effect_solver.h"
#include "gelex/algo/reml/operators.h"
#include "gelex/algo/reml/reml_buffer.h"
#include "gelex/algo/reml/summary.h"
#include "gelex/freq/model.h"

namespace gelex
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

auto wald_p_onesided(double z) noexcept -> double
{
    return 0.5 * std::erfc(z / std::numbers::sqrt2);
}

auto wald_p_twosided(double z) noexcept -> double
{
    return std::erfc(std::abs(z) / std::numbers::sqrt2);
}

auto compute_variance_se(FreqState& state, const RemlBuffer& buffer) -> void
{
    // se(σ) = sqrt(diag(-H⁻¹))
    // variance component order: residual, random[0..]
    Eigen::VectorXd se = (-buffer.hess_inv.diagonal()).array().sqrt();

    // A component clamped to the constraint floor sits on the parameter-space
    // boundary, where the AI Hessian gives no valid Wald SE; report NaN.
    const double floor = constraint_floor(buffer.phenotype_variance());
    constexpr double nan = std::numeric_limits<double>::quiet_NaN();
    auto at_boundary
        = [floor](double variance) noexcept { return variance <= floor; };

    Eigen::Index idx = 0;

    // residual
    state.residual().variance_se
        = at_boundary(state.residual().variance) ? nan : se(idx);
    ++idx;

    // random effects
    for (auto& r : state.random())
    {
        r.variance_se = at_boundary(r.variance) ? nan : se(idx);
        ++idx;
    }
}

auto compute_variance_ratio(FreqState& state, const RemlBuffer& buffer) -> void
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
    auto n_comp = buffer.hess_inv.rows();

    const double floor = constraint_floor(buffer.phenotype_variance());
    constexpr double nan = std::numeric_limits<double>::quiet_NaN();

    for (size_t ri = 0; ri < state.random().size(); ++ri)
    {
        auto& r = state.random()[ri];
        Eigen::Index r_idx = 1 + static_cast<Eigen::Index>(ri);

        r.variance_ratio = r.variance / sum_var;

        // A boundary component has no valid Wald SE (see compute_variance_se);
        // keep the ratio point estimate but drop its delta-method SE.
        if (r.variance <= floor)
        {
            r.variance_ratio_se = nan;
            continue;
        }

        // gradient for delta method:
        // ∂ratio/∂σ_i = -σ_r / (Σσ)² for i ≠ r
        // ∂ratio/∂σ_r = (Σσ - σ_r) / (Σσ)² for i == r
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

        double var_ratio = grad.dot(-buffer.hess_inv * grad);
        r.variance_ratio_se = std::sqrt(std::max(0.0, var_ratio));
    }
}

auto assemble_reml_fit(
    const FreqModel& model,
    FreqState& state,
    RemlBuffer& buffer,
    double loglike,
    bool converged,
    size_t iter_count) -> RemlFit
{
    compute_fixed_effects(model, state, buffer);
    compute_random_effects(model, state, buffer);
    compute_variance_se(state, buffer);
    compute_variance_ratio(state, buffer);

    std::vector<VarianceComponent> components;
    components.reserve(state.random().size());
    for (size_t i = 0; i < state.random().size(); ++i)
    {
        const auto& r = state.random()[i];
        components.push_back(
            {.name = model.random()[i].name,
             .variance = r.variance,
             .variance_se = r.variance_se,
             .variance_ratio = r.variance_ratio,
             .variance_ratio_se = r.variance_ratio_se});
    }

    RemlSummary summary{
        .loglike = loglike,
        .converged = converged,
        .iter_count = iter_count,
        .random = std::move(components),
        .residual_variance = state.residual().variance,
        .residual_variance_se = state.residual().variance_se};

    // Materialize P = V^{-1} - ViX * XtViX_inv * ViX' in buffer.V's memory
    // so downstream GWAS can use a single dense GEMM per SNP chunk.
    buffer.V.noalias()
        -= buffer.ViX * buffer.XtViX_inv * buffer.ViX.transpose();

    return RemlFit{
        .summary = std::move(summary),
        .operators = GwasOperators{
            .P = std::move(buffer.V),
            .Py = std::move(buffer.Py),
            .Vp = state.Vp()}};
}

}  // namespace gelex
