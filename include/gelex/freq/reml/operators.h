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
