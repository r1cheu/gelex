// Copyright 2026 RuLei Chen
// SPDX-License-Identifier: Apache-2.0

#include "gelex/freq/reml/effect_solver.h"

#include <cstddef>

#include "gelex/freq/model.h"
#include "gelex/freq/reml/reml_buffer.h"

namespace gelex
{

auto compute_fixed_effects(
    const FreqModel& model,
    FreqState& state,
    const RemlBuffer& buffer) -> void
{
    // β = inv_XtViX * ViX' * y = (X'V⁻¹X)⁻¹ * X' * V⁻¹ * y
    // se(β) = sqrt(diag(inv_XtViX))
    state.fixed().coeffs.noalias()
        = buffer.XtViX_inv * (buffer.ViX.transpose() * model.phenotype());
    state.fixed().se = buffer.XtViX_inv.diagonal().array().sqrt();
}

auto compute_random_effects(
    const FreqModel& model,
    FreqState& state,
    const RemlBuffer& buffer) -> void
{
    for (size_t i = 0; i < model.random().size(); ++i)
    {
        const auto& effect = model.random()[i];
        auto& effect_state = state.random()[i];

        effect_state.blup.noalias()
            = effect.K * buffer.Py * effect_state.variance;
    }
}

}  // namespace gelex
