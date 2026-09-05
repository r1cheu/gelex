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

#ifndef GELEX_BAYES_GENETIC_DETAIL_COEFFICIENT_LIKELIHOOD_H_
#define GELEX_BAYES_GENETIC_DETAIL_COEFFICIENT_LIKELIHOOD_H_

#include <Eigen/Core>

#include "gelex/bayes/basic_state.h"
#include "gelex/bayes/genotype/projection.h"
#include "gelex/bayes/stats/quadratic_log_kernel.h"

namespace gelex::detail
{

[[nodiscard]] inline auto make_coefficient_likelihood(
    const bayes::GeneticProjection& projection,
    Eigen::Index marker,
    double current_coefficient,
    const ResidualState& residual) -> QuadraticLogKernel
{
    const double quadratic = projection.xtx_diag()(marker);
    const double linear = projection.dot(marker, residual.adjusted_response)
                          + (quadratic * current_coefficient);

    return gelex::make_coefficient_likelihood(
        quadratic, linear, residual.variance);
}

}  // namespace gelex::detail

#endif  // GELEX_BAYES_GENETIC_DETAIL_COEFFICIENT_LIKELIHOOD_H_
