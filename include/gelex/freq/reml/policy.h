// Copyright 2026 RuLei Chen
// SPDX-License-Identifier: Apache-2.0

#ifndef GELEX_FREQ_REML_POLICY_H_
#define GELEX_FREQ_REML_POLICY_H_

#include <Eigen/Core>

#include "gelex/freq/model.h"
#include "gelex/freq/reml/reml_buffer.h"

namespace gelex
{

struct EMPolicy
{
    static auto apply(
        const gelex::FreqModel& model,
        const gelex::FreqState& state,
        RemlBuffer& buffer) -> Eigen::VectorXd;
};

struct AIPolicy
{
    // AI-REML search direction delta = -H^{-1} * grad, evaluated at the point
    // whose V^{-1}/P/Py are currently held in buffer. The caller owns the
    // step length (see the Armijo backtracking in Estimator::fit).
    static auto direction(const gelex::FreqModel& model, RemlBuffer& buffer)
        -> Eigen::VectorXd;
};

}  // namespace gelex

#endif  // GELEX_FREQ_REML_POLICY_H_
