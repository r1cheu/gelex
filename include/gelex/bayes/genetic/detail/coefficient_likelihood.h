// Copyright 2026 RuLei Chen
// SPDX-License-Identifier: Apache-2.0

#ifndef GELEX_BAYES_GENETIC_DETAIL_COEFFICIENT_LIKELIHOOD_H_
#define GELEX_BAYES_GENETIC_DETAIL_COEFFICIENT_LIKELIHOOD_H_

#include <Eigen/Core>

#include "gelex/bayes/genotype/projection.h"
#include "gelex/bayes/state.h"
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
