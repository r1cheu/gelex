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

#include "gelex/algo/infer/vi/steps/fixed.h"

#include <Eigen/Core>

namespace gelex::vi
{

auto FixedStep::step() -> void
{
    const auto& X = deps_.effect.X;
    const auto& XtX_diag = deps_.effect.XtX_diag;
    auto& coeffs = deps_.state.coeffs;
    auto& y_adj = deps_.residual.y_adj;

    for (Eigen::Index i = 0; i < coeffs.size(); ++i)
    {
        const double old_i = coeffs(i);
        const auto col = X.col(i);
        const double norm = XtX_diag(i);

        const double rhs = col.dot(y_adj) + (norm * old_i);
        const double post_mean = rhs / norm;

        coeffs(i) = post_mean;
        y_adj.array() += (old_i - post_mean) * col.array();
    }
}

}  // namespace gelex::vi
