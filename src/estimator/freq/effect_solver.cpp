#include "gelex/estimator/freq/effect_solver.h"

#include <Eigen/Cholesky>

#include "gelex/model/freq/model_new.h"
#include "gelex/optim/optimizer_state.h"

namespace gelex::effect_solver
{

auto compute_fixed_effects(
    const FreqModel& model,
    FreqState& state,
    const OptimizerState& opt_state) -> void
{
    const auto& x = model.fixed().X;
    const auto& y = model.phenotype();

    // (X'V⁻¹X)⁻¹
    // note: opt_state.v contains V⁻¹ after v_inv_logdet
    Eigen::LLT<Eigen::MatrixXd> llt(opt_state.tx_vinv_x);
    Eigen::MatrixXd tx_vinv_x_inv = llt.solve(
        Eigen::MatrixXd::Identity(
            opt_state.tx_vinv_x.rows(), opt_state.tx_vinv_x.cols()));

    // β = (X'V⁻¹X)⁻¹ * X' * V⁻¹ * y
    state.fixed().coeff = tx_vinv_x_inv * (x.transpose() * opt_state.v * y);

    // se(β) = sqrt(diag((X'V⁻¹X)⁻¹))
    state.fixed().se = tx_vinv_x_inv.diagonal().array().sqrt();
}

auto compute_random_effects(
    const FreqModel& model,
    FreqState& state,
    const OptimizerState& opt_state) -> void
{
    // random effects: blup = K * Py * σ
    for (size_t i = 0; i < model.random().size(); ++i)
    {
        const auto& effect = model.random()[i];
        auto& effect_state = state.random()[i];

        effect_state.blup.noalias()
            = effect.K * opt_state.proj_y * effect_state.variance;
    }

    // genetic effects: ebv = K * Py * σ
    for (size_t i = 0; i < model.genetic().size(); ++i)
    {
        const auto& effect = model.genetic()[i];
        auto& effect_state = state.genetic()[i];

        effect_state.ebv.noalias()
            = effect.K * opt_state.proj_y * effect_state.variance;
    }
}

}  // namespace gelex::effect_solver
