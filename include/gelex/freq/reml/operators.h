// Copyright 2026 RuLei Chen
// SPDX-License-Identifier: Apache-2.0

#ifndef GELEX_FREQ_REML_OPERATORS_H_
#define GELEX_FREQ_REML_OPERATORS_H_

#include <Eigen/Core>

namespace gelex
{

// REML projection operators consumed by the association scan. P is the mixed
// -model projection matrix, materialized (n x n) by reusing the V^{-1} buffer
// at the end of Estimator::fit so per-chunk scans reduce to a single dense
// GEMM; Py is P applied to the phenotype, Vp the phenotypic variance.
struct GwasOperators
{
    Eigen::MatrixXd P;
    Eigen::VectorXd Py;
    double Vp{};

    [[nodiscard]] auto n_samples() const -> Eigen::Index { return P.rows(); }
};

}  // namespace gelex

#endif  // GELEX_FREQ_REML_OPERATORS_H_
