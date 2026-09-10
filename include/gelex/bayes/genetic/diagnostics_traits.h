// Copyright 2026 RuLei Chen
// SPDX-License-Identifier: Apache-2.0

#ifndef GELEX_BAYES_GENETIC_DIAGNOSTICS_TRAITS_H_
#define GELEX_BAYES_GENETIC_DIAGNOSTICS_TRAITS_H_

#include <Eigen/Core>
#include <cstddef>
#include <ranges>
#include <span>
#include <string>
#include <string_view>
#include <type_traits>
#include <utility>
#include <variant>
#include <vector>

#include "gelex/bayes/draws_diagnostics.h"
#include "gelex/bayes/genetic/types.h"
#include "gelex/bayes/stats/diagnostics.h"
#include "gelex/io/csc_reader.h"
#include "gelex/io/dense_reader.h"
#include "gelex/namespace.h"

GELEX_NAMESPACE_BEGIN(gelex)

// Diagnostic factories read dense payloads from `dense` and marker-level
// sparse payloads from `sparse`; the mirror of DrawWriters.
struct DrawReaders
{
    const DenseReader& dense;
    const CscReader& sparse;
};

template <VarianceLayout Kind>
using marker_variance_diagnostics_t = std::conditional_t<
    Kind == VarianceLayout::Pooled,
    ChainDiagnostics,
    std::monostate>;

template <VarianceLayout Kind>
auto diagnose_marker_variance(
    DrawReaders readers,
    std::string_view identifier,
    double prob) -> marker_variance_diagnostics_t<Kind>
{
    if constexpr (Kind == VarianceLayout::Pooled)
    {
        return detail::diagnose_scalar(readers.dense, identifier, prob);
    }
    else
    {
        return {};
    }
}

template <MixtureWeightUpdate Update>
using probability_diagnostics_t = std::conditional_t<
    Update == MixtureWeightUpdate::Enabled,
    ChainDiagnostics,
    std::monostate>;

template <MixtureWeightUpdate Update>
using probabilities_diagnostics_t = std::conditional_t<
    Update == MixtureWeightUpdate::Enabled,
    std::vector<ChainDiagnostics>,
    std::monostate>;

template <MixtureWeightUpdate Update>
auto diagnose_probability(
    DrawReaders readers,
    std::string_view identifier,
    double prob) -> probability_diagnostics_t<Update>
{
    if constexpr (Update == MixtureWeightUpdate::Enabled)
    {
        return detail::diagnose_scalar(readers.dense, identifier, prob);
    }
    else
    {
        return {};
    }
}

template <MixtureWeightUpdate Update>
auto diagnose_probabilities(
    DrawReaders readers,
    std::string_view identifier,
    double prob) -> probabilities_diagnostics_t<Update>
{
    if constexpr (Update == MixtureWeightUpdate::Enabled)
    {
        return detail::diagnose_rows(readers.dense, identifier, prob);
    }
    else
    {
        return {};
    }
}

// One diagnosed parameter, addressed by its payload id (or a derived-quantity
// id) and, for multi-row payloads, its row.
struct DiagnosticEntry
{
    std::string id;
    Eigen::Index index;
    ChainDiagnostics stats;
};

inline auto append_entry(
    std::vector<DiagnosticEntry>& out,
    std::string id,
    const ChainDiagnostics& stats) -> void
{
    out.push_back({.id = std::move(id), .index = 0, .stats = stats});
}

inline auto append_entries(
    std::vector<DiagnosticEntry>& out,
    std::string id,
    std::span<const ChainDiagnostics> stats) -> void
{
    for (const auto& [index, value] : std::views::enumerate(stats))
    {
        out.push_back({.id = id, .index = index, .stats = value});
    }
}

// Absent payloads (Unpooled variance, Disabled mixture weights) yield nothing.
inline auto append_entry(
    std::vector<DiagnosticEntry>& /*out*/,
    std::string /*id*/,
    std::monostate /*stats*/) -> void
{
}

inline auto append_entries(
    std::vector<DiagnosticEntry>& /*out*/,
    std::string /*id*/,
    std::monostate /*stats*/) -> void
{
}

// Marker coefficient draws as (n_markers, n_draws); both maps alias the
// reader's file and must not outlive it.
template <CoefficientLayout Layout>
auto read_coefficients(DrawReaders readers, std::string_view identifier)
{
    if constexpr (Layout == CoefficientLayout::Dense)
    {
        return readers.dense.to_map<double>(identifier);
    }
    else
    {
        return readers.sparse.to_map<double>(identifier);
    }
}

GELEX_NAMESPACE_END(gelex)

#endif  // GELEX_BAYES_GENETIC_DIAGNOSTICS_TRAITS_H_
