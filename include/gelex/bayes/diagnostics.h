// Copyright 2026 RuLei Chen
// SPDX-License-Identifier: Apache-2.0

#ifndef GELEX_BAYES_DIAGNOSTICS_H_
#define GELEX_BAYES_DIAGNOSTICS_H_

#include <Eigen/Core>
#include <cstddef>
#include <fmt/format.h>
#include <span>
#include <string>
#include <string_view>
#include <type_traits>
#include <utility>
#include <vector>

#include "gelex/bayes/draws.h"
#include "gelex/bayes/draws_diagnostics.h"
#include "gelex/bayes/genetic/diagnostics_traits.h"
#include "gelex/bayes/genetic/factory.h"
#include "gelex/bayes/genotype/gebv.h"
#include "gelex/bayes/mode_values.h"
#include "gelex/bayes/model.h"
#include "gelex/bayes/serialization_ids.h"
#include "gelex/bayes/stats/diagnostics.h"
#include "gelex/bayes/variance/heritability.h"
#include "gelex/genetic_mode.h"
#include "gelex/infra/var.h"
#include "gelex/io/csc_reader.h"
#include "gelex/io/dense_reader.h"
#include "gelex/namespace.h"

GELEX_NAMESPACE_BEGIN(gelex)

// Post-run diagnostics of every model term written by BayesDraws.
template <typename GeneticPrior>
class BayesDiagnostics
{
   public:
    static constexpr GeneticModeSet modes = GeneticPrior::modes;
    using genetic_draws_type = genetic_draws_t<genetic_state_t<GeneticPrior>>;
    using genetic_diagnostics_type = genetic_diagnostics_t<genetic_draws_type>;
    using genetic_variance_type
        = HomogeneousModeValues<modes, GeneticVarianceDiagnostics>;

    BayesDiagnostics(
        std::vector<ChainDiagnostics> fixed,
        std::vector<RandomEffectDiagnostics> random,
        genetic_diagnostics_type genetic,
        genetic_variance_type genetic_variance,
        GeneticVarianceDiagnostics total_genetic_variance,
        ChainDiagnostics residual)
        : fixed_{std::move(fixed)},
          random_{std::move(random)},
          genetic_{std::move(genetic)},
          genetic_variance_{std::move(genetic_variance)},
          total_genetic_variance_{total_genetic_variance},
          residual_{residual}
    {
    }

    [[nodiscard]] auto fixed() const noexcept
        -> const std::vector<ChainDiagnostics>&
    {
        return fixed_;
    }
    [[nodiscard]] auto random() const noexcept
        -> const std::vector<RandomEffectDiagnostics>&
    {
        return random_;
    }
    [[nodiscard]] auto genetic() const noexcept
        -> const genetic_diagnostics_type&
    {
        return genetic_;
    }
    template <GeneticMode Mode>
    [[nodiscard]] auto genetic_variance() const noexcept
        -> const GeneticVarianceDiagnostics&
    {
        return genetic_variance_.template get<Mode>();
    }
    [[nodiscard]] auto total_genetic_variance() const noexcept
        -> const GeneticVarianceDiagnostics&
    {
        return total_genetic_variance_;
    }
    [[nodiscard]] auto residual() const noexcept -> const ChainDiagnostics&
    {
        return residual_;
    }

   private:
    std::vector<ChainDiagnostics> fixed_;
    std::vector<RandomEffectDiagnostics> random_;
    genetic_diagnostics_type genetic_;
    genetic_variance_type genetic_variance_;
    GeneticVarianceDiagnostics total_genetic_variance_;
    ChainDiagnostics residual_;
};

GELEX_NAMESPACE_BEGIN(detail)

template <GeneticModeSet Modes>
struct ExplainedVarianceDraws
{
    HomogeneousModeValues<Modes, Eigen::RowVectorXd> mode;
    Eigen::RowVectorXd total;
};

// Population variance of every draw's GEBV per mode and of their sum, each
// (1, n_draws). One draw's GEBV lives only in a per-thread scratch vector.
template <typename GeneticDraws>
auto explained_variance_draws(
    DrawReaders readers,
    const bayes::GeneticDesign& design)
    -> ExplainedVarianceDraws<GeneticDraws::modes>
{
    constexpr auto modes = GeneticDraws::modes;
    const auto coefficients = generate_mode_values<modes>(
        [&]<GeneticMode Mode>()
        {
            using draws_type =
                typename GeneticDraws::template mode_value_type<Mode>;
            return read_coefficients<draws_type::coefficient_layout>(
                readers,
                fmt::format("{}/{}", genetic_id<Mode>, coefficients_id));
        });
    const Eigen::Index n_draws
        = coefficients.template get<modes.at(0)>().cols();
    const Eigen::Index n_individuals = design.rows();
    ExplainedVarianceDraws<modes> result{
        .mode
        = generate_mode_values<modes>([&]<GeneticMode /*Mode*/>()
                                      { return Eigen::RowVectorXd(n_draws); }),
        .total = Eigen::RowVectorXd(n_draws)};

#pragma omp parallel default(none) \
    shared(coefficients, design, result, n_draws, n_individuals)
    {
        auto gebv = generate_mode_values<modes>(
            [&]<GeneticMode /*Mode*/>()
            { return Eigen::VectorXd(n_individuals); });
#pragma omp for
        for (Eigen::Index draw = 0; draw < n_draws; ++draw)
        {
            gebv.for_each(
                [&]<GeneticMode Mode>(Eigen::VectorXd& target)
                {
                    gebv_draw(
                        design.projection(Mode),
                        coefficients.template get<Mode>(),
                        draw,
                        target);
                    result.mode.template get<Mode>()(draw)
                        = vecvar(target, VarNormType::Population);
                });
            if constexpr (modes.size() == 1)
            {
                result.total(draw)
                    = result.mode.template get<modes.at(0)>()(draw);
            }
            else
            {
                result.total(draw) = vecvar(
                    gebv.template get<GeneticMode::A>()
                        + gebv.template get<GeneticMode::D>(),
                    VarNormType::Population);
            }
        }
    }
    return result;
}

GELEX_NAMESPACE_END(detail)

template <typename GeneticPrior>
[[nodiscard]] auto read_diagnostics(
    std::string_view draws_path,
    const BayesModel& model,
    double prob = 0.95) -> BayesDiagnostics<GeneticPrior>
{
    using diagnostics_type = BayesDiagnostics<GeneticPrior>;
    using genetic_draws_type = typename diagnostics_type::genetic_draws_type;
    const DenseReader dense{draws_path};
    const CscReader sparse{sparse_draws_path(draws_path)};
    const DrawReaders readers{.dense = dense, .sparse = sparse};

    std::vector<RandomEffectDiagnostics> random;
    random.reserve(model.random().size());
    for (const auto& design : model.random())
    {
        random.push_back(diagnose_random(dense, design.name(), prob));
    }

    const auto explained = detail::explained_variance_draws<genetic_draws_type>(
        readers, model.genetic());
    const auto residual_payload = dense.to_map<double>(residual_variance_id);
    const Eigen::RowVectorXd residual_variance = residual_payload.row(0);

    const auto diagnose_variance
        = [&](const Eigen::RowVectorXd& variance) -> GeneticVarianceDiagnostics
    {
        return {
            .explained_variance = diagnose_chain(variance, prob),
            .heritability = diagnose_chain(
                heritability_draws(
                    variance, explained.total, residual_variance),
                prob)};
    };
    auto genetic_variance = transform_mode_values(
        explained.mode,
        [&]<GeneticMode /*Mode*/>(const Eigen::RowVectorXd& variance)
        { return diagnose_variance(variance); });

    return diagnostics_type{
        diagnose_fixed(dense, prob),
        std::move(random),
        make_diagnostics(
            std::type_identity<genetic_draws_type>{}, readers, prob),
        std::move(genetic_variance),
        diagnose_variance(explained.total),
        diagnose_residual(dense, prob)};
}

// Every diagnosed parameter in a fixed order: fixed coefficients, each random
// effect's coefficients then variance, the genetic family payloads per mode,
// each mode's explained variance and heritability, the A + D totals when both
// modes are present, and the residual variance.
template <typename GeneticPrior>
[[nodiscard]] auto diagnostic_entries(
    const BayesDiagnostics<GeneticPrior>& diagnostics)
    -> std::vector<DiagnosticEntry>
{
    constexpr auto modes = BayesDiagnostics<GeneticPrior>::modes;
    std::vector<DiagnosticEntry> entries;
    append_entries(
        entries, std::string{fixed_coefficients_id}, diagnostics.fixed());
    for (const auto& random : diagnostics.random())
    {
        append_entries(
            entries, random_coefficients_id(random.name), random.coefficients);
        append_entry(entries, random_variance_id(random.name), random.variance);
    }
    append_diagnostic_entries(entries, diagnostics.genetic());

    const auto append_variance =
        [&](std::string_view prefix, const GeneticVarianceDiagnostics& variance)
    {
        append_entry(
            entries,
            fmt::format("{}/{}", prefix, explained_variance_id),
            variance.explained_variance);
        append_entry(
            entries,
            fmt::format("{}/{}", prefix, heritability_id),
            variance.heritability);
    };
    [&]<std::size_t... Index>(std::index_sequence<Index...>)
    {
        (append_variance(
             genetic_id<modes.at(Index)>,
             diagnostics.template genetic_variance<modes.at(Index)>()),
         ...);
    }(std::make_index_sequence<modes.size()>{});
    if constexpr (modes.size() > 1)
    {
        append_variance(total_genetic_id, diagnostics.total_genetic_variance());
    }
    append_entry(
        entries, std::string{residual_variance_id}, diagnostics.residual());
    return entries;
}

// Tab-separated table with one row per entry: id, index, mean, sd, median,
// hpdi_lower, hpdi_upper, ess, mcse, split_rhat.
auto write_diagnostics(
    std::string_view path,
    std::span<const DiagnosticEntry> entries) -> void;

GELEX_NAMESPACE_END(gelex)

#endif  // GELEX_BAYES_DIAGNOSTICS_H_
