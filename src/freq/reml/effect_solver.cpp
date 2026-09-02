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
