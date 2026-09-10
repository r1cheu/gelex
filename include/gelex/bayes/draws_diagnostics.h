// Copyright 2026 RuLei Chen
// SPDX-License-Identifier: Apache-2.0

#ifndef GELEX_BAYES_DRAWS_DIAGNOSTICS_H_
#define GELEX_BAYES_DRAWS_DIAGNOSTICS_H_

#include <string_view>
#include <vector>

#include "gelex/bayes/stats/diagnostics.h"
#include "gelex/io/dense_reader.h"

namespace gelex
{

namespace detail
{

// One diagnostic per payload row (each row is a parameter's chain).
auto diagnose_rows(
    const DenseReader& draws,
    std::string_view identifier,
    double prob) -> std::vector<ChainDiagnostics>;

// Diagnostic of a single-row payload; any other row count throws.
auto diagnose_scalar(
    const DenseReader& draws,
    std::string_view identifier,
    double prob) -> ChainDiagnostics;

}  // namespace detail

struct RandomEffectDiagnostics
{
    // One entry per level, in payload row order.
    std::vector<ChainDiagnostics> coefficients;
    ChainDiagnostics variance;
};

// Diagnostics of every fixed-effect coefficient, in payload row order.
auto diagnose_fixed(const DenseReader& draws, double prob = 0.95)
    -> std::vector<ChainDiagnostics>;

// Diagnostics of the random effect `name`: its level coefficients and
// variance component.
auto diagnose_random(
    const DenseReader& draws,
    std::string_view name,
    double prob = 0.95) -> RandomEffectDiagnostics;

// Diagnostics of the residual variance.
auto diagnose_residual(const DenseReader& draws, double prob = 0.95)
    -> ChainDiagnostics;

}  // namespace gelex

#endif  // GELEX_BAYES_DRAWS_DIAGNOSTICS_H_
