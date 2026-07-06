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

#include "gelex/algo/reml/variance_calculator.h"

#include <Eigen/Core>
#include <cmath>
#include <numbers>
#include <ranges>
#include <stdexcept>

#include <Eigen/src/misc/lapacke.h>
#include <Eigen/Cholesky>

#include "gelex/algo/reml/optimizer_state.h"
#include "gelex/freq/model.h"

namespace gelex
{

auto compute_v(
    const gelex::FreqModel& model,
    const gelex::FreqState& state,
    Eigen::Ref<Eigen::MatrixXd> v) -> void
{
    v.setZero();

    // residual: I * sigma_e
    v.diagonal().array() += state.residual().variance;

    auto compute_v = [&](auto& effects, auto& states)
    {
        for (auto&& [effect, state] : std::views::zip(effects, states))
        {
            v += effect.K * state.variance;
        }
    };
    compute_v(model.random(), state.random());
}

auto v_inv_logdet(Eigen::Ref<Eigen::MatrixXd> v) -> double
{
    auto n = static_cast<lapack_int>(v.cols());
    lapack_int info = 0;

    // 1. Cholesky decomposition: V = L * L^T
    info = LAPACKE_dpotrf(LAPACK_COL_MAJOR, 'L', n, v.data(), n);
    if (info != 0)
    {
        throw std::runtime_error("V matrix is not positive definite");
    }

    // 2. Compute log|V| = 2 * sum(log(diag(L)))
    // GCC false positive: deep inlining of Eigen expression templates
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wmaybe-uninitialized"
    double logdet = 2.0 * v.diagonal().array().log().sum();
#pragma GCC diagnostic pop

    // 3. Compute inverse from Cholesky factor: V^{-1}
    info = LAPACKE_dpotri(LAPACK_COL_MAJOR, 'L', n, v.data(), n);
    if (info != 0)
    {
        throw std::runtime_error("Error computing inverse of V");
    }

    // 4. Symmetrize (potri only fills lower triangle)
    v.triangularView<Eigen::StrictlyUpper>()
        = v.triangularView<Eigen::StrictlyLower>().transpose();

    return logdet;
}

auto compute_proj(const gelex::FreqModel& model, OptimizerState& state) -> void
{
    const auto& x = model.fixed().X;
    const auto& y = model.phenotype();

    // ViX = V^{-1} X (v holds V^{-1} after v_inv_logdet)
    state.ViX.noalias() = state.V * x;

    // XtViX = X' V^{-1} X = X' ViX (local only; stored form is its inverse)
    Eigen::MatrixXd XtViX = x.transpose() * state.ViX;

    Eigen::LLT<Eigen::MatrixXd> llt(XtViX);
    if (llt.info() != Eigen::Success)
    {
        throw std::runtime_error("X'V^{-1}X is not positive definite");
    }

    // log|X'V^{-1}X| = 2 * sum(log(diag(L)))
    Eigen::MatrixXd L = llt.matrixL();
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wmaybe-uninitialized"
    state.logdet_xvx = 2.0 * L.diagonal().array().log().sum();
#pragma GCC diagnostic pop

    // inv_XtViX = (X' V^{-1} X)^{-1}
    state.XtViX_inv
        = llt.solve(Eigen::MatrixXd::Identity(XtViX.rows(), XtViX.cols()));

    // Py = V^{-1} y - ViX inv_XtViX (ViX' y)
    state.Py.noalias() = state.V * y;
    state.Py.noalias()
        -= state.ViX * (state.XtViX_inv * (state.ViX.transpose() * y));
}

auto compute_loglike(const gelex::FreqModel& model, const OptimizerState& state)
    -> double
{
    // logL = -0.5 * ((n-p)*log(2π) + log|V| + log|X'V^{-1}X| + y'Py)
    double ypy = model.phenotype().dot(state.Py);
    auto dof
        = static_cast<double>(model.num_individuals() - model.fixed().X.cols());
    double constant = dof * std::log(2.0 * std::numbers::pi);
    return -0.5 * (constant + state.logdet_v + state.logdet_xvx + ypy);
}

}  // namespace gelex
