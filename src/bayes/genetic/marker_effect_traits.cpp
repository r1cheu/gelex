// Copyright 2026 RuLei Chen
// SPDX-License-Identifier: Apache-2.0

#include "gelex/bayes/genetic/marker_effect_traits.h"

#include <Eigen/Core>
#include <Eigen/SparseCore>
#include <algorithm>
#include <fmt/format.h>
#include <string>
#include <string_view>
#include <utility>

#include "gelex/exception.h"
#include "gelex/genetic_mode.h"
#include "gelex/infra/var.h"
#include "gelex/io/csc_reader.h"

namespace gelex
{

namespace detail
{

namespace
{

auto require_draws(Eigen::Index draw_count) -> void
{
    if (draw_count < 2)
    {
        throw GelexException(
            fmt::format(
                "marker effects: need at least 2 draws, got {}", draw_count));
    }
}

}  // namespace

auto summarize_coefficients(
    const Eigen::Ref<const Eigen::MatrixXd>& coefficients)
    -> MarkerCoefficientSummary
{
    require_draws(coefficients.cols());
    return {
        .mean = coefficients.rowwise().mean(),
        .sd = matvar<1>(coefficients, VarNormType::Sample).cwiseSqrt()};
}

auto summarize_coefficients(
    const CscReader::sparse_map_type<double>& coefficients)
    -> MarkerCoefficientSummary
{
    const Eigen::Index draw_count = coefficients.cols();
    require_draws(draw_count);
    Eigen::VectorXd sum = Eigen::VectorXd::Zero(coefficients.rows());
    Eigen::VectorXd sum_sq = Eigen::VectorXd::Zero(coefficients.rows());
    for (Eigen::Index draw = 0; draw < draw_count; ++draw)
    {
        for (CscReader::sparse_map_type<double>::InnerIterator it(
                 coefficients, draw);
             it;
             ++it)
        {
            const auto row = static_cast<Eigen::Index>(it.row());
            sum(row) += it.value();
            sum_sq(row) += it.value() * it.value();
        }
    }
    const auto n = static_cast<double>(draw_count);
    Eigen::VectorXd mean = sum / n;
    Eigen::VectorXd sd
        = ((sum_sq.array() - n * mean.array().square()) / (n - 1.0))
              .max(0.0)
              .sqrt();
    return {.mean = std::move(mean), .sd = std::move(sd)};
}

}  // namespace detail

MarkerEffectTable::MarkerEffectTable(Eigen::Index rows) : rows_{rows}
{
    if (rows < 0)
    {
        throw GelexException(
            "MarkerEffectTable: row count must not be "
            "negative");
    }
}

auto MarkerEffectTable::add(std::string name, Eigen::VectorXd values) -> void
{
    if (values.size() != rows_)
    {
        throw GelexException(
            fmt::format(
                "MarkerEffectTable: column \"{}\" has {} rows, expected {}",
                name,
                values.size(),
                rows_));
    }
    if (std::ranges::find(names_, name) != names_.end())
    {
        throw GelexException(
            fmt::format("MarkerEffectTable: duplicate column \"{}\"", name));
    }
    names_.push_back(std::move(name));
    columns_.push_back(std::move(values));
}

auto MarkerEffectTable::column(std::string_view name) const
    -> const Eigen::VectorXd&
{
    const auto it = std::ranges::find(names_, name);
    if (it == names_.end())
    {
        throw GelexException(
            fmt::format("MarkerEffectTable: no column \"{}\"", name));
    }
    return columns_[static_cast<std::size_t>(it - names_.begin())];
}

MarkerPveScale::MarkerPveScale(
    const Eigen::Ref<const Eigen::RowVectorXd>& col_var,
    double denominator)
    : scale_{col_var.transpose() / denominator}
{
    if (!(denominator > 0.0))
    {
        throw GelexException(
            fmt::format(
                "MarkerPveScale: denominator must be positive, got {}",
                denominator));
    }
}

auto MarkerPveScale::pve(const Eigen::Ref<const Eigen::VectorXd>& mean) const
    -> Eigen::VectorXd
{
    if (mean.size() != scale_.size())
    {
        throw GelexException(
            fmt::format(
                "MarkerPveScale: {} coefficients but {} marker variances",
                mean.size(),
                scale_.size()));
    }
    return mean.array().square() * scale_.array();
}

auto append_marker_columns(
    MarkerEffectTable& out,
    const CoefficientMarkerEffects& effects,
    GeneticMode mode,
    const MarkerPveScale& pve) -> void
{
    out.add(fmt::format("BETA_{}", mode), effects.coefficients.mean);
    out.add(fmt::format("SE_{}", mode), effects.coefficients.sd);
    out.add(fmt::format("PVE_{}", mode), pve.pve(effects.coefficients.mean));
}

auto append_marker_columns(
    MarkerEffectTable& out,
    const MixtureMarkerEffects& effects,
    GeneticMode mode,
    const MarkerPveScale& pve) -> void
{
    append_marker_columns(
        out, CoefficientMarkerEffects{effects.coefficients}, mode, pve);
    out.add(fmt::format("PIP_{}", mode), effects.pip);
}

}  // namespace gelex
