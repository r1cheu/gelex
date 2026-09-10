// Copyright 2026 RuLei Chen
// SPDX-License-Identifier: Apache-2.0

#ifndef GELEX_BAYES_VARIANCE_HERITABILITY_H_
#define GELEX_BAYES_VARIANCE_HERITABILITY_H_

#include <Eigen/Core>

#include "gelex/bayes/stats/diagnostics.h"

namespace gelex
{

// Draw-wise explained / (total genetic + residual); every input is
// (1, n_draws).
auto heritability_draws(
    const Eigen::Ref<const Eigen::RowVectorXd>& explained,
    const Eigen::Ref<const Eigen::RowVectorXd>& total_genetic,
    const Eigen::Ref<const Eigen::RowVectorXd>& residual) -> Eigen::RowVectorXd;

struct GeneticVarianceDiagnostics
{
    ChainDiagnostics explained_variance;
    ChainDiagnostics heritability;
};

}  // namespace gelex

#endif  // GELEX_BAYES_VARIANCE_HERITABILITY_H_
