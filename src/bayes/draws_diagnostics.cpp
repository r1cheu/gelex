// Copyright 2026 RuLei Chen
// SPDX-License-Identifier: Apache-2.0

#include "gelex/bayes/draws_diagnostics.h"

#include <Eigen/Core>
#include <cstddef>
#include <fmt/format.h>
#include <string>
#include <string_view>
#include <vector>

#include "gelex/bayes/serialization_ids.h"
#include "gelex/bayes/stats/diagnostics.h"
#include "gelex/exception.h"
#include "gelex/io/dense_reader.h"

namespace gelex
{

namespace detail
{

auto diagnose_rows(
    const DenseReader& draws,
    std::string_view identifier,
    double prob) -> std::vector<ChainDiagnostics>
{
    const auto payload = draws.to_map<double>(identifier);
    std::vector<ChainDiagnostics> result;
    result.reserve(static_cast<std::size_t>(payload.rows()));
    for (Eigen::Index row = 0; row < payload.rows(); ++row)
    {
        result.push_back(diagnose_chain(payload.row(row), prob));
    }
    return result;
}

auto diagnose_scalar(
    const DenseReader& draws,
    std::string_view identifier,
    double prob) -> ChainDiagnostics
{
    const auto payload = draws.to_map<double>(identifier);
    if (payload.rows() != 1)
    {
        throw GelexException(
            fmt::format(
                "payload \"{}\" must hold one row, got {}",
                identifier,
                payload.rows()));
    }
    return diagnose_chain(payload.row(0), prob);
}

}  // namespace detail

auto diagnose_fixed(const DenseReader& draws, double prob)
    -> std::vector<ChainDiagnostics>
{
    return detail::diagnose_rows(draws, fixed_coefficients_id, prob);
}

auto diagnose_random(
    const DenseReader& draws,
    std::string_view name,
    double prob) -> RandomEffectDiagnostics
{
    return RandomEffectDiagnostics{
        .name = std::string{name},
        .coefficients
        = detail::diagnose_rows(draws, random_coefficients_id(name), prob),
        .variance
        = detail::diagnose_scalar(draws, random_variance_id(name), prob),
    };
}

auto diagnose_residual(const DenseReader& draws, double prob)
    -> ChainDiagnostics
{
    return detail::diagnose_scalar(draws, residual_variance_id, prob);
}

}  // namespace gelex
