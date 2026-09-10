// Copyright 2026 RuLei Chen
// SPDX-License-Identifier: Apache-2.0

#ifndef GELEX_BAYES_GENETIC_MARKER_EFFECT_TRAITS_H_
#define GELEX_BAYES_GENETIC_MARKER_EFFECT_TRAITS_H_

#include <Eigen/Core>
#include <Eigen/SparseCore>
#include <concepts>
#include <cstdint>
#include <span>
#include <string>
#include <string_view>
#include <vector>

#include "gelex/bayes/genetic/diagnostics_traits.h"
#include "gelex/bayes/genetic/types.h"
#include "gelex/genetic_mode.h"
#include "gelex/io/csc_reader.h"
#include "gelex/namespace.h"

GELEX_NAMESPACE_BEGIN(gelex)

// Posterior mean and sample standard deviation of every marker coefficient.
struct MarkerCoefficientSummary
{
    Eigen::VectorXd mean;
    Eigen::VectorXd sd;
};

GELEX_NAMESPACE_BEGIN(detail)

auto summarize_coefficients(
    const Eigen::Ref<const Eigen::MatrixXd>& coefficients)
    -> MarkerCoefficientSummary;

auto summarize_coefficients(
    const CscReader::sparse_map_type<double>& coefficients)
    -> MarkerCoefficientSummary;

GELEX_NAMESPACE_END(detail)

template <CoefficientLayout Layout>
auto summarize_coefficients(DrawReaders readers, std::string_view identifier)
    -> MarkerCoefficientSummary
{
    return detail::summarize_coefficients(
        read_coefficients<Layout>(readers, identifier));
}

// Fraction of draws whose stored class satisfies `active`, per marker. CSC
// payloads omit zero classes, so `active(0)` must be false.
template <std::predicate<std::uint8_t> Active>
auto inclusion_probability(
    const CscReader& sparse,
    std::string_view identifier,
    Active active) -> Eigen::VectorXd
{
    const auto assignments = sparse.to_map<std::uint8_t>(identifier);
    Eigen::VectorXd counts = Eigen::VectorXd::Zero(assignments.rows());
    for (Eigen::Index draw = 0; draw < assignments.cols(); ++draw)
    {
        for (CscReader::sparse_map_type<std::uint8_t>::InnerIterator it(
                 assignments, draw);
             it;
             ++it)
        {
            if (active(it.value()))
            {
                counts(static_cast<Eigen::Index>(it.row())) += 1.0;
            }
        }
    }
    return counts / static_cast<double>(assignments.cols());
}

// Ordered, uniquely named numeric columns of equal length.
class MarkerEffectTable
{
   public:
    explicit MarkerEffectTable(Eigen::Index rows);

    auto add(std::string name, Eigen::VectorXd values) -> void;

    [[nodiscard]] auto rows() const noexcept -> Eigen::Index { return rows_; }
    [[nodiscard]] auto names() const noexcept -> std::span<const std::string>
    {
        return names_;
    }
    [[nodiscard]] auto columns() const noexcept
        -> std::span<const Eigen::VectorXd>
    {
        return columns_;
    }
    [[nodiscard]] auto column(std::string_view name) const
        -> const Eigen::VectorXd&;

   private:
    Eigen::Index rows_;
    std::vector<std::string> names_;
    std::vector<Eigen::VectorXd> columns_;
};

// Turns one mode's posterior mean coefficients into per-marker explained
// variance fractions: var(x_j) * mean_j^2 / denominator.
class MarkerPveScale
{
   public:
    MarkerPveScale(
        const Eigen::Ref<const Eigen::RowVectorXd>& col_var,
        double denominator);

    [[nodiscard]] auto pve(const Eigen::Ref<const Eigen::VectorXd>& mean) const
        -> Eigen::VectorXd;

   private:
    Eigen::VectorXd scale_;
};

struct CoefficientMarkerEffects
{
    MarkerCoefficientSummary coefficients;
};

struct MixtureMarkerEffects
{
    MarkerCoefficientSummary coefficients;
    Eigen::VectorXd pip;
};

// BETA_<mode>, SE_<mode>, PVE_<mode>.
auto append_marker_columns(
    MarkerEffectTable& out,
    const CoefficientMarkerEffects& effects,
    GeneticMode mode,
    const MarkerPveScale& pve) -> void;

// BETA_<mode>, SE_<mode>, PVE_<mode>, PIP_<mode>.
auto append_marker_columns(
    MarkerEffectTable& out,
    const MixtureMarkerEffects& effects,
    GeneticMode mode,
    const MarkerPveScale& pve) -> void;

GELEX_NAMESPACE_END(gelex)

#endif  // GELEX_BAYES_GENETIC_MARKER_EFFECT_TRAITS_H_
