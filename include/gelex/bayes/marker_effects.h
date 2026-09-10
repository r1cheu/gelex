// Copyright 2026 RuLei Chen
// SPDX-License-Identifier: Apache-2.0

#ifndef GELEX_BAYES_MARKER_EFFECTS_H_
#define GELEX_BAYES_MARKER_EFFECTS_H_

#include <Eigen/Core>
#include <string_view>
#include <type_traits>

#include "gelex/bayes/draws.h"
#include "gelex/bayes/genetic/diagnostics_traits.h"
#include "gelex/bayes/genetic/factory.h"
#include "gelex/bayes/genetic/marker_effect_traits.h"
#include "gelex/bayes/genotype/design.h"
#include "gelex/bayes/genotype/gebv.h"
#include "gelex/bayes/mode_values.h"
#include "gelex/bayes/model.h"
#include "gelex/bayes/serialization_ids.h"
#include "gelex/genetic_mode.h"
#include "gelex/infra/var.h"
#include "gelex/io/csc_reader.h"
#include "gelex/io/dense_reader.h"
#include "gelex/namespace.h"

GELEX_NAMESPACE_BEGIN(gelex)

GELEX_NAMESPACE_BEGIN(detail)

template <typename Effects>
[[nodiscard]] auto mode_marker_effects(const Effects& effects) -> decltype(auto)
{
    if constexpr (requires { effects.mode_values(); })
    {
        return (effects.mode_values());
    }
    else
    {
        return (effects);
    }
}

GELEX_NAMESPACE_END(detail)

// Per-marker posterior summaries of every genetic mode, scaled by the
// posterior-mean genetic plus residual variance: BETA/SE/PVE per mode, PIP
// where the family tracks inclusion, and the A + D total PVE when both modes
// are present.
template <typename GeneticPrior>
[[nodiscard]] auto read_marker_effects(
    std::string_view draws_path,
    const BayesModel& model) -> MarkerEffectTable
{
    constexpr auto modes = GeneticPrior::modes;
    using genetic_draws_type = genetic_draws_t<genetic_state_t<GeneticPrior>>;
    const DenseReader dense{draws_path};
    const CscReader sparse{sparse_draws_path(draws_path)};
    const DrawReaders readers{.dense = dense, .sparse = sparse};
    const auto& design = model.genetic();

    const auto effects = make_marker_effects(
        std::type_identity<genetic_draws_type>{}, readers);
    const auto& mode_effects = detail::mode_marker_effects(effects);

    Eigen::VectorXd gebv = Eigen::VectorXd::Zero(design.rows());
    Eigen::VectorXd scratch(design.rows());
    mode_effects.for_each(
        [&]<GeneticMode Mode>(const auto& mode_effect)
        {
            gebv_draw(
                design.projection(Mode),
                mode_effect.coefficients.mean,
                0,
                scratch);
            gebv += scratch;
        });
    const double residual_variance
        = dense.to_map<double>(residual_variance_id).row(0).mean();
    const double denominator
        = vecvar(gebv, VarNormType::Population) + residual_variance;

    const auto scales = generate_mode_values<modes>(
        [&]<GeneticMode Mode>()
        {
            return MarkerPveScale{
                design.projection(Mode).col_var(), denominator};
        });

    MarkerEffectTable table{design.cols()};
    append_marker_columns(table, effects, scales);
    if constexpr (modes.size() > 1)
    {
        const auto& additive = design.projection(GeneticMode::A);
        const auto& dominance = design.projection(GeneticMode::D);
        const Eigen::VectorXd beta_a
            = mode_effects.template get<GeneticMode::A>().coefficients.mean;
        const Eigen::VectorXd beta_d
            = mode_effects.template get<GeneticMode::D>().coefficients.mean;
        const Eigen::VectorXd covariance
            = additive.col_covariance(dominance).transpose();
        table.add(
            "PVE",
            (beta_a.array().square() * additive.col_var().transpose().array()
             + beta_d.array().square() * dominance.col_var().transpose().array()
             + 2.0 * beta_a.array() * beta_d.array() * covariance.array())
                / denominator);
    }
    return table;
}

// Tab-separated table: CHR, SNP, BP, A1, A2 and A1FREQ from the design's
// marker metadata, followed by every column of `table`.
auto write_marker_effects(
    std::string_view path,
    const bayes::GeneticDesign& design,
    const MarkerEffectTable& table) -> void;

GELEX_NAMESPACE_END(gelex)

#endif  // GELEX_BAYES_MARKER_EFFECTS_H_
