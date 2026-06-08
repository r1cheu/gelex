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

#include "gelex/algo/reml/effect_solver.h"
#include <cstddef>

#include "gelex/algo/reml/optimizer_state.h"
#include "gelex/freq/model.h"

namespace gelex::reml
{

auto compute_fixed_effects(
    const FreqModel& model,
    FreqState& state,
    const OptimizerState& opt_state) -> void
{
    // β = inv_XtViX * ViX' * y = (X'V⁻¹X)⁻¹ * X' * V⁻¹ * y
    // se(β) = sqrt(diag(inv_XtViX))
    state.fixed().coeffs.noalias()
        = opt_state.XtViX_inv * (opt_state.ViX.transpose() * model.phenotype());
    state.fixed().se = opt_state.XtViX_inv.diagonal().array().sqrt();
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
            = effect.K * opt_state.Py * effect_state.variance;
    }

    // genetic effects: ebv = K * Py * σ
    for (size_t i = 0; i < model.genetic().size(); ++i)
    {
        const auto& effect = model.genetic()[i];
        auto& effect_state = state.genetic()[i];

        effect_state.ebv.noalias()
            = effect.K * opt_state.Py * effect_state.variance;
    }
}

}  // namespace gelex::reml
