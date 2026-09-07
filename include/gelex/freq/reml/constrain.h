// Copyright 2026 RuLei Chen
// SPDX-License-Identifier: Apache-2.0

#ifndef GELEX_FREQ_REML_CONSTRAIN_H_
#define GELEX_FREQ_REML_CONSTRAIN_H_

#include <Eigen/Dense>

namespace gelex
{
constexpr double constraint_scale = 1e-6;

// Lower limit that variance components are clamped to. A component at or below
// this floor is statistically indistinguishable from zero, so its Wald SE has
// no meaning.
inline auto constraint_floor(double y_variance) noexcept -> double
{
    return y_variance * constraint_scale;
}

// Clamps negative variance components to a small positive limit, redistributing
// the borrowed mass across the unconstrained ones. Returns a boundary mask
// flagging which components were clamped, so the caller can judge reliability
// (mask.count()) and mark boundary components whose Wald test is invalid.
auto constrain(Eigen::Ref<Eigen::VectorXd> varcmp, double y_variance)
    -> Eigen::ArrayX<bool>;
}  // namespace gelex

#endif  // GELEX_FREQ_REML_CONSTRAIN_H_
