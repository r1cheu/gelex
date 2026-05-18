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

#include "gelex/algo/infer/vi/steps/random.h"

#include <cstddef>

#include <fmt/format.h>
#include <Eigen/Core>

#include "gelex/exception.h"
#include "gelex/infra/stats/conjugate_prior.h"

namespace gelex::vi
{

auto RandomStep::make(const Context& ctx) -> RandomStep
{
    const auto& effects = ctx.model.random();
    const auto& specs = ctx.method.randoms;
    auto& states = ctx.state.random();
    if (effects.size() != specs.size() || effects.size() != states.size())
    {
        throw GelexException(
            fmt::format(
                "random block size mismatch: model={}, method.random={}, "
                "state={}",
                effects.size(),
                specs.size(),
                states.size()));
    }
    return RandomStep{Deps{
        .effects = std::span{effects},
        .specs = std::span{specs},
        .states = std::span{states},
        .residual = ctx.state.residual(),
    }};
}

auto RandomStep::step() -> void
{
    auto& y_adj = deps_.residual.y_adj;
    const double residual_variance = deps_.residual.variance;

    for (std::size_t block = 0; block < deps_.effects.size(); ++block)
    {
        const auto& effect = deps_.effects[block];
        const auto& spec = deps_.specs[block];
        auto& state = deps_.states[block];

        auto& coeffs = state.coeffs;
        const auto& X = effect.X;
        const auto& XtX_diag = effect.XtX_diag;
        const double sigma = state.variance;

        const Eigen::VectorXd inv_scaler
            = 1.0 / (XtX_diag.array() + residual_variance / sigma);

        for (Eigen::Index i = 0; i < coeffs.size(); ++i)
        {
            const double old_i = coeffs(i);
            const auto col = X.col(i);
            const double norm = XtX_diag(i);

            const double rhs = col.dot(y_adj) + (norm * old_i);
            const double post_mean = rhs * inv_scaler(i);

            coeffs(i) = post_mean;
            y_adj.array() += (old_i - post_mean) * col.array();
        }

        gelex::stats::detail::ScaledInvChiSq chi_squared{
            {.nu = spec.prior.degrees_of_freedom(), .s2 = spec.prior.scale()}};
        chi_squared.compute(coeffs.squaredNorm(), coeffs.size());
        state.variance = chi_squared.expected_value();
    }
}

}  // namespace gelex::vi
