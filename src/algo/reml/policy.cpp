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

#include "gelex/algo/reml/policy.h"

#include <Eigen/Core>
#include <Eigen/Dense>

#include "gelex/algo/reml/reml_buffer.h"
#include "gelex/algo/reml/variance_calculator.h"
#include "gelex/freq/model.h"

namespace gelex
{

auto EMPolicy::apply(
    const gelex::FreqModel& model,
    const gelex::FreqState& state,
    RemlBuffer& buffer) -> Eigen::VectorXd
{
    Eigen::VectorXd sigma = collect_variance_components(state);
    Eigen::VectorXd sigma_sq = sigma.array().square();
    auto n = static_cast<double>(buffer.num_individuals());

    // residual: K = I
    // sigma_new = (sigma^2 * Py'Py - sigma^2 * tr(P) + sigma * n) / n
    double py_py = buffer.Py.squaredNorm();
    double tr_p = buffer.trace_proj();
    sigma(0) = (sigma_sq(0) * py_py - sigma_sq(0) * tr_p + sigma(0) * n) / n;

    Eigen::Index idx = 1;

    auto compute_sigma = [&](const auto& effects)
    {
        for (const auto& effect : effects)
        {
            double py_k_py = buffer.Py.dot(effect.K * buffer.Py);
            double tr_pk = buffer.trace_proj_k(effect.K);
            sigma(idx) = (sigma_sq(idx) * py_k_py - sigma_sq(idx) * tr_pk
                          + sigma(idx) * n)
                         / n;

            ++idx;
        }
    };
    compute_sigma(model.random());
    return sigma;
}

auto AIPolicy::direction(const gelex::FreqModel& model, RemlBuffer& buffer)
    -> Eigen::VectorXd
{
    auto n_comp = static_cast<Eigen::Index>(1 + model.random().size());

    // 1. Compute dvpy: dV/dsigma_i * P * y for each component
    // residual: dV/dsigma_0 = I, so dvpy(:,0) = Py
    buffer.dvpy.col(0) = buffer.Py;

    Eigen::Index idx = 1;

    auto compute_dvpy = [&](const auto& effects)
    {
        for (const auto& effect : effects)
        {
            buffer.dvpy.col(idx++).noalias() = effect.K * buffer.Py;
        }
    };
    compute_dvpy(model.random());

    // 2. Compute first gradient
    // grad(i) = -0.5 * (tr(P * dV/dsigma_i) - Py' * dV/dsigma_i * Py)
    // = -0.5 * (tr(P * K_i) - Py' * K_i * Py)
    // residual: K_0 = I
    buffer.first_grad(0)
        = -0.5 * (buffer.trace_proj() - buffer.Py.dot(buffer.dvpy.col(0)));

    idx = 1;
    auto compute_first_grad = [&](const auto& effects)
    {
        for (const auto& effect : effects)
        {
            double tr_pk = buffer.trace_proj_k(effect.K);
            double py_k_py = buffer.Py.dot(buffer.dvpy.col(idx));
            buffer.first_grad(idx) = -0.5 * (tr_pk - py_k_py);
            ++idx;
        }
    };
    compute_first_grad(model.random());

    // 3. Compute AI Hessian: H(i,j) = -0.5 * dvpy(:,i)' * P * dvpy(:,j)
    // P * dvpy = V^{-1} * dvpy - ViX * inv_XtViX * (ViX' * dvpy)
    Eigen::MatrixXd p_dvpy = buffer.V * buffer.dvpy;
    Eigen::MatrixXd vix_dvpy = buffer.ViX.transpose() * buffer.dvpy;
    p_dvpy.noalias() -= buffer.ViX * (buffer.XtViX_inv * vix_dvpy);

    // Only compute upper triangle since Hessian is symmetric
    Eigen::MatrixXd hess(n_comp, n_comp);
    for (Eigen::Index i = 0; i < n_comp; ++i)
    {
        for (Eigen::Index j = i; j < n_comp; ++j)
        {
            hess(i, j) = -0.5 * buffer.dvpy.col(i).dot(p_dvpy.col(j));
            if (i != j)
            {
                hess(j, i) = hess(i, j);
            }
        }
    }

    // 4. Compute direction: delta = -H^{-1} * grad
    buffer.hess_inv = hess.completeOrthogonalDecomposition().pseudoInverse();
    return -buffer.hess_inv * buffer.first_grad;
}

}  // namespace gelex
