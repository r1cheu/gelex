// Copyright 2026 RuLei Chen
// SPDX-License-Identifier: Apache-2.0

#include "gelex/freq/design.h"

#include <Eigen/Core>

#include "gelex/data/fixed_design.h"

namespace gelex::freq
{

FixedState::FixedState(const gelex::FixedDesign& design)
    : coeffs(Eigen::VectorXd::Zero(design.X().cols())),
      se(Eigen::VectorXd::Zero(design.X().cols()))
{
}

RandomState::RandomState(const RandomDesign& design)
    : blup(Eigen::VectorXd::Zero(design.K.rows()))
{
}

}  // namespace gelex::freq
