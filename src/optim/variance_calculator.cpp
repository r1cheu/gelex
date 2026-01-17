#include "gelex/optim/variance_calculator.h"

#include <ranges>
#include <stdexcept>

#include <Eigen/src/misc/lapacke.h>
#include <Eigen/Cholesky>

#include "gelex/model/freq/model.h"
#include "gelex/optim/optimizer_state.h"

namespace gelex::variance_calculator
{

auto compute_v(
    const FreqModel& model,
    const FreqState& state,
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
    compute_v(model.genetic(), state.genetic());
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
    double logdet = 2.0 * v.diagonal().array().log().sum();

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

auto compute_proj(const FreqModel& model, OptimizerState& state) -> void
{
    const auto& x = model.fixed().X;

    // v now contains V^{-1}
    // vinv_x = V^{-1} * X
    Eigen::MatrixXd vinv_x = state.v * x;

    // tx_vinv_x = X' * V^{-1} * X
    state.tx_vinv_x.noalias() = x.transpose() * vinv_x;

    // solve (X'V^{-1}X)^{-1} * (V^{-1}X)'
    Eigen::LLT<Eigen::MatrixXd> llt_xvx(state.tx_vinv_x);
    if (llt_xvx.info() != Eigen::Success)
    {
        throw std::runtime_error("X'V^{-1}X is not positive definite");
    }

    // P = V^{-1} - V^{-1}*X * (X'V^{-1}X)^{-1} * X'*V^{-1}
    // P = V^{-1} - vinv_x * (tx_vinv_x)^{-1} * vinv_x'
    Eigen::MatrixXd xvx_inv_xv = llt_xvx.solve(vinv_x.transpose());
    state.proj.noalias() = state.v - vinv_x * xvx_inv_xv;

    // Py = P * y
    state.proj_y.noalias() = state.proj * model.phenotype();
}

auto compute_loglike(const FreqModel& model, const OptimizerState& state)
    -> double
{
    // logL = -0.5 * (log|V| + log|X'V^{-1}X| + y'Py)

    // log|X'V^{-1}X|
    Eigen::LLT<Eigen::MatrixXd> llt_xvx(state.tx_vinv_x);
    Eigen::MatrixXd L_xvx = llt_xvx.matrixL();
    double logdet_xvx = 2.0 * L_xvx.diagonal().array().log().sum();

    // y'Py
    double ypy = model.phenotype().dot(state.proj_y);

    return -0.5 * (state.logdet_v + logdet_xvx + ypy);
}

}  // namespace gelex::variance_calculator
