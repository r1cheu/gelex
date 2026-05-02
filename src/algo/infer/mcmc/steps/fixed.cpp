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

#include "gelex/algo/infer/mcmc/steps/fixed.h"

#include <Eigen/Core>

#include "gelex/infra/stats/conjugate_prior.h"

namespace gelex::mcmc
{

auto FixedStep::step() -> void
{
    normal_.reset();

    const auto& X = deps_.effect.X;
    const auto& XtX_diag = deps_.effect.XtX_diag;
    auto& coeffs = deps_.state.coeffs;
    auto& y_adj = deps_.residual.y_adj;
    const double residual_variance = deps_.residual.variance;

    // Flat prior on fixed effects: posterior is the limit prior_var -> inf,
    // so Params are computed in closed form and fed to NormalSampler::draw().
    for (Eigen::Index i = 0; i < coeffs.size(); ++i)
    {
        const double old_i = coeffs(i);
        const auto col = X.col(i);
        const double norm = XtX_diag(i);

        const double rhs = col.dot(y_adj) + (norm * old_i);
        const stats::NormalSampler<double>::Params posterior{
            .mean = rhs / norm,
            .var = residual_variance / norm,
        };
        const double new_i = normal_.draw(posterior, deps_.rng);
        coeffs(i) = new_i;

        y_adj.array() += (old_i - new_i) * col.array();
    }
}

}  // namespace gelex::mcmc
