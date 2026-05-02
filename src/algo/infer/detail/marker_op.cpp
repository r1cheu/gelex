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

#include "gelex/algo/infer/detail/marker_op.h"

#include <cmath>
#include <cstdint>
#include <limits>
#include <span>

#include <Eigen/Core>

namespace gelex::infer::detail
{

// NOLINTBEGIN(cppcoreguidelines-pro-bounds-pointer-arithmetic)
__attribute__((target_clones("default,avx2,avx512f"))) auto
apply_marker_update_impl(
    double* y_adj,  // NOLINT(bugprone-easily-swappable-parameters)
    double* u,
    std::span<Eigen::VectorXd> component_u,
    const double* col,
    Eigen::Index n,
    const MarkerTransition& tx) -> void
{
    const double diff = tx.old_value - tx.new_value;
    const bool diff_zero
        = std::abs(diff) <= std::numeric_limits<double>::epsilon();
    const int8_t oc = tx.old_class;
    const int8_t nc = tx.new_class;
    const bool no_mix = component_u.empty() || (oc <= 0 && nc <= 0);

    if (diff_zero && (no_mix || oc == nc))
    {
        return;
    }

    const double* __restrict c = col;
    double* __restrict ya = y_adj;
    double* __restrict uu = u;

    if (no_mix)
    {
        for (Eigen::Index j = 0; j < n; ++j)
        {
            const double cj = c[j];
            ya[j] += diff * cj;
            uu[j] -= diff * cj;
        }
        return;
    }

    if (oc <= 0)
    {
        double* __restrict cu = component_u[nc - 1].data();
        const double nv = tx.new_value;
        for (Eigen::Index j = 0; j < n; ++j)
        {
            const double cj = c[j];
            ya[j] += diff * cj;
            uu[j] -= diff * cj;
            cu[j] += nv * cj;
        }
        return;
    }

    if (nc <= 0)
    {
        double* __restrict cu = component_u[oc - 1].data();
        const double ov = tx.old_value;
        for (Eigen::Index j = 0; j < n; ++j)
        {
            const double cj = c[j];
            ya[j] += diff * cj;
            uu[j] -= diff * cj;
            cu[j] -= ov * cj;
        }
        return;
    }

    if (oc == nc)
    {
        // same component, value changed: component_u[k-1] += (-diff) * col
        double* __restrict cu = component_u[oc - 1].data();
        for (Eigen::Index j = 0; j < n; ++j)
        {
            const double cj = c[j];
            ya[j] += diff * cj;
            uu[j] -= diff * cj;
            cu[j] -= diff * cj;
        }
        return;
    }

    double* __restrict cu_o = component_u[oc - 1].data();
    double* __restrict cu_n = component_u[nc - 1].data();
    const double ov = tx.old_value;
    const double nv = tx.new_value;
    for (Eigen::Index j = 0; j < n; ++j)
    {
        const double cj = c[j];
        ya[j] += diff * cj;
        uu[j] -= diff * cj;
        cu_o[j] -= ov * cj;
        cu_n[j] += nv * cj;
    }
}
// NOLINTEND(cppcoreguidelines-pro-bounds-pointer-arithmetic)

}  // namespace gelex::infer::detail
