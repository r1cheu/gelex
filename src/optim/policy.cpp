#include "gelex/optim/policy.h"

#include <Eigen/Dense>

#include "gelex/model/freq/model.h"
#include "gelex/optim/optimizer.h"
#include "gelex/optim/optimizer_state.h"

namespace gelex
{

auto EMPolicy::apply(
    const FreqModel& model,
    const FreqState& state,
    OptimizerState& opt_state) -> Eigen::VectorXd
{
    Eigen::VectorXd sigma = collect_variance_components(state);
    Eigen::VectorXd sigma_sq = sigma.array().square();
    auto n = static_cast<double>(opt_state.num_individuals());

    // residual: K = I
    // sigma_new = (sigma^2 * Py'Py - sigma^2 * tr(P) + sigma * n) / n
    double py_py = opt_state.proj_y.squaredNorm();
    double tr_p = opt_state.proj.trace();
    sigma(0) = (sigma_sq(0) * py_py - sigma_sq(0) * tr_p + sigma(0) * n) / n;

    Eigen::Index idx = 1;

    auto compute_sigma = [&](const auto& effects)
    {
        for (const auto& effect : effects)
        {
            double py_k_py = opt_state.proj_y.dot(effect.K * opt_state.proj_y);
            double tr_pk = (opt_state.proj * effect.K).trace();
            sigma(idx) = (sigma_sq(idx) * py_k_py - sigma_sq(idx) * tr_pk
                          + sigma(idx) * n)
                         / n;

            ++idx;
        }
    };
    compute_sigma(model.random());
    compute_sigma(model.genetic());
    return sigma;
}

auto AIPolicy::apply(
    const FreqModel& model,
    const FreqState& state,
    OptimizerState& opt_state) -> Eigen::VectorXd
{
    Eigen::VectorXd sigma = collect_variance_components(state);
    auto n_comp = sigma.size();
    auto n = opt_state.num_individuals();

    // 1. Compute dvpy: dV/dsigma_i * P * y for each component
    // residual: dV/dsigma_0 = I, so dvpy(:,0) = Py
    opt_state.dvpy.col(0) = opt_state.proj_y;

    Eigen::Index idx = 1;

    auto compute_dvpy = [&](const auto& effects)
    {
        for (const auto& effect : effects)
        {
            opt_state.dvpy.col(idx++).noalias() = effect.K * opt_state.proj_y;
        }
    };
    compute_dvpy(model.random());
    compute_dvpy(model.genetic());

    // 2. Compute first gradient
    // grad(i) = -0.5 * (tr(P * dV/dsigma_i) - Py' * dV/dsigma_i * Py)
    //         = -0.5 * (tr(P * K_i) - Py' * K_i * Py)
    // residual: K_0 = I
    opt_state.first_grad(0) = -0.5
                              * (opt_state.proj.trace()
                                 - opt_state.proj_y.dot(opt_state.dvpy.col(0)));

    idx = 1;
    auto compute_first_grad = [&](const auto& effects)
    {
        for (const auto& effect : effects)
        {
            double tr_pk = (opt_state.proj * effect.K).trace();
            double py_k_py = opt_state.proj_y.dot(opt_state.dvpy.col(idx));
            opt_state.first_grad(idx) = -0.5 * (tr_pk - py_k_py);
            ++idx;
        }
    };
    compute_first_grad(model.random());
    compute_first_grad(model.genetic());

    // 3. Compute AI Hessian: H(i,j) = -0.5 * dvpy(:,i)' * P * dvpy(:,j)
    // Precompute P * dvpy to avoid redundant matrix-vector products
    Eigen::MatrixXd p_dvpy = opt_state.proj * opt_state.dvpy;

    // Only compute upper triangle since Hessian is symmetric
    Eigen::MatrixXd hess(n_comp, n_comp);
    for (Eigen::Index i = 0; i < n_comp; ++i)
    {
        for (Eigen::Index j = i; j < n_comp; ++j)
        {
            hess(i, j) = -0.5 * opt_state.dvpy.col(i).dot(p_dvpy.col(j));
            if (i != j)
            {
                hess(j, i) = hess(i, j);
            }
        }
    }

    // 4. Compute update: delta = -H^{-1} * grad
    opt_state.hess_inv = hess.completeOrthogonalDecomposition().pseudoInverse();
    Eigen::VectorXd delta = -opt_state.hess_inv * opt_state.first_grad;

    return sigma + delta;
}

}  // namespace gelex
