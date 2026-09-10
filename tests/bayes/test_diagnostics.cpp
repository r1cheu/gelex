// Copyright 2026 RuLei Chen
// SPDX-License-Identifier: Apache-2.0

#include <Eigen/Core>
#include <algorithm>
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <cstddef>
#include <fstream>
#include <string>
#include <type_traits>
#include <vector>

#include "gelex/bayes/builtin_method.h"
#include "gelex/bayes/diagnostics.h"
#include "gelex/bayes/draws.h"
#include "gelex/bayes/draws_diagnostics.h"
#include "gelex/bayes/genotype/gebv.h"
#include "gelex/bayes/prior.h"
#include "gelex/bayes/state.h"
#include "gelex/bayes/variance/budget.h"
#include "gelex/bayes/variance/heritability.h"
#include "gelex/genetic_mode.h"
#include "gelex/infra/var.h"
#include "gelex/io/dense_reader.h"

#include "bayes_model_fixture.h"
#include "file_fixture.h"

namespace
{

constexpr auto mode_ad = gelex::GeneticMode::A | gelex::GeneticMode::D;
constexpr std::size_t draw_count = 4;

using Catch::Matchers::WithinAbs;

}  // namespace

TEST_CASE(
    "Bayes diagnostics cover every term and combine A + D variance",
    "[bayes][diagnostics]")
{
    gelex::test::FileFixture fixture;
    const auto path = (fixture.get_test_dir() / "rr.draws").string();
    const auto model = gelex::test::make_random_effect_model(mode_ad);
    const auto prior = gelex::make_prior(
        gelex::BuiltinBayesRecipe<mode_ad, gelex::BayesMethod::RR>{
            gelex::VarianceBudget{
                {.additive = 0.4, .dominance = 0.1, .random = 0.1}}},
        model);
    using prior_type = std::remove_cvref_t<decltype(prior)>::genetic_prior_type;
    auto state = gelex::make_state(prior, model);
    {
        auto draws = gelex::BayesDraws{state, model, path, draw_count};
        for (std::size_t draw = 0; draw < draw_count; ++draw)
        {
            const double value = static_cast<double>(draw + 1);
            state.fixed().coefficients = Eigen::VectorXd{{value}};
            state.random()[0].coefficients = Eigen::VectorXd{{value, -value}};
            state.random()[0].variance = value;
            state.residual().variance = 2.0 * value;
            state.genetic().get<gelex::GeneticMode::A>().transition(
                0, 0.5 * value);
            state.genetic().get<gelex::GeneticMode::A>().transition(1, -value);
            state.genetic().get<gelex::GeneticMode::D>().transition(
                1, 0.25 * value);
            draws << state;
        }
        draws.close();
    }

    const auto result = gelex::read_diagnostics<prior_type>(path, model, 0.9);
    const gelex::DenseReader reader{path};

    SECTION("fixed, random and residual match the term functions")
    {
        REQUIRE(result.fixed().size() == 1);
        REQUIRE(
            result.fixed()[0].mean == gelex::diagnose_fixed(reader)[0].mean);
        REQUIRE(result.random().size() == 1);
        REQUIRE(result.random()[0].coefficients.size() == 2);
        REQUIRE(result.random()[0].variance.mean == 2.5);
        REQUIRE(result.residual().mean == 5.0);
    }

    SECTION("genetic diagnostics carry the pooled variance per mode")
    {
        const auto& additive
            = result.genetic().get<gelex::GeneticMode::A>().variance;
        REQUIRE(
            additive.mean
            == gelex::diagnose_chain(
                   reader.to_map<double>("genetic/A/variance").row(0), 0.9)
                   .mean);
    }

    SECTION("total genetic variance is var(gebv_A + gebv_D)")
    {
        const auto additive = reader.to_map<double>("genetic/A/coefficients");
        const auto dominance = reader.to_map<double>("genetic/D/coefficients");
        Eigen::RowVectorXd explained_a(draw_count);
        Eigen::RowVectorXd explained_d(draw_count);
        Eigen::RowVectorXd total(draw_count);
        Eigen::VectorXd gebv_a(model.num_individuals());
        Eigen::VectorXd gebv_d(model.num_individuals());
        for (Eigen::Index draw = 0; draw < additive.cols(); ++draw)
        {
            gelex::gebv_draw(
                model.genetic().projection(gelex::GeneticMode::A),
                additive,
                draw,
                gebv_a);
            gelex::gebv_draw(
                model.genetic().projection(gelex::GeneticMode::D),
                dominance,
                draw,
                gebv_d);
            explained_a(draw)
                = gelex::vecvar(gebv_a, gelex::VarNormType::Population);
            explained_d(draw)
                = gelex::vecvar(gebv_d, gelex::VarNormType::Population);
            total(draw) = gelex::vecvar(
                gebv_a + gebv_d, gelex::VarNormType::Population);
        }

        REQUIRE_THAT(
            result.genetic_variance<gelex::GeneticMode::A>()
                .explained_variance.mean,
            WithinAbs(explained_a.mean(), 1e-12));
        REQUIRE_THAT(
            result.genetic_variance<gelex::GeneticMode::D>()
                .explained_variance.mean,
            WithinAbs(explained_d.mean(), 1e-12));
        REQUIRE_THAT(
            result.total_genetic_variance().explained_variance.mean,
            WithinAbs(total.mean(), 1e-12));
        REQUIRE(
            result.total_genetic_variance().explained_variance.mean
            != explained_a.mean() + explained_d.mean());

        const Eigen::RowVectorXd residual
            = reader.to_map<double>("residual/variance").row(0);
        const Eigen::RowVectorXd heritability
            = gelex::heritability_draws(total, total, residual);
        REQUIRE_THAT(
            result.total_genetic_variance().heritability.mean,
            WithinAbs(heritability.mean(), 1e-12));
        REQUIRE(result.total_genetic_variance().heritability.mean > 0.0);
        REQUIRE(result.total_genetic_variance().heritability.mean < 1.0);
    }

    SECTION("entries list every parameter under its payload id")
    {
        const auto entries = gelex::diagnostic_entries(result);
        std::vector<std::string> ids;
        for (const auto& entry : entries)
        {
            ids.push_back(entry.id);
        }
        const std::string random_name{model.random()[0].name()};
        const std::vector<std::string> expected{
            "fixed/coefficients",
            "random/" + random_name + "/coefficients",
            "random/" + random_name + "/coefficients",
            "random/" + random_name + "/variance",
            "genetic/A/variance",
            "genetic/D/variance",
            "genetic/A/explained_variance",
            "genetic/A/heritability",
            "genetic/D/explained_variance",
            "genetic/D/heritability",
            "genetic/total/explained_variance",
            "genetic/total/heritability",
            "residual/variance"};
        REQUIRE(ids == expected);
        REQUIRE(entries[2].index == 1);
        REQUIRE(
            entries[2].stats.mean == result.random()[0].coefficients[1].mean);
        REQUIRE(entries.back().stats.mean == result.residual().mean);
        REQUIRE(
            entries[11].stats.mean
            == result.total_genetic_variance().heritability.mean);
    }

    SECTION("write_diagnostics emits a header and one row per entry")
    {
        const auto entries = gelex::diagnostic_entries(result);
        const auto summary_path
            = (fixture.get_test_dir() / "rr.summary").string();
        gelex::write_diagnostics(summary_path, entries);

        std::ifstream file{summary_path};
        std::vector<std::string> lines;
        for (std::string line; std::getline(file, line);)
        {
            lines.push_back(line);
        }
        REQUIRE(lines.size() == entries.size() + 1);
        REQUIRE(
            lines[0]
            == "id\tindex\tmean\tsd\tmedian\thpdi_lower\thpdi_upper\tess\tmcse"
               "\tsplit_rhat");
        REQUIRE(lines[1].starts_with("fixed/coefficients\t0\t2.5\t"));
        REQUIRE(lines.back().starts_with("residual/variance\t0\t5\t"));
    }
}

TEST_CASE(
    "Bayes diagnostics dispatch the joint genetic family",
    "[bayes][diagnostics]")
{
    gelex::test::FileFixture fixture;
    const auto path = (fixture.get_test_dir() / "cd.draws").string();
    const auto model = gelex::test::make_random_effect_model(mode_ad);
    const auto prior = gelex::make_prior(
        gelex::BuiltinBayesRecipe<mode_ad, gelex::BayesMethod::CD>{
            gelex::VarianceBudget{
                {.additive = 0.4, .dominance = 0.1, .random = 0.1}}},
        model);
    using prior_type = std::remove_cvref_t<decltype(prior)>::genetic_prior_type;
    auto state = gelex::make_state(prior, model);
    {
        auto draws = gelex::BayesDraws{state, model, path, draw_count};
        for (std::size_t draw = 0; draw < draw_count; ++draw)
        {
            const double value = static_cast<double>(draw + 1);
            state.residual().variance = value;
            state.genetic().get<gelex::GeneticMode::A>().transition(0, value);
            state.genetic().get<gelex::GeneticMode::D>().transition(
                1, -0.5 * value);
            state.genetic().joint().transition(0, 1);
            state.genetic().joint().transition(1, 2);
            draws << state;
        }
        draws.close();
    }

    const auto result = gelex::read_diagnostics<prior_type>(path, model);
    const auto entries = gelex::diagnostic_entries(result);
    REQUIRE(
        std::ranges::count(
            entries, "genetic/joint/probabilities", &gelex::DiagnosticEntry::id)
        == 4);
    REQUIRE(
        std::ranges::count(
            entries,
            "genetic/D/annotation_coefficients",
            &gelex::DiagnosticEntry::id)
        == 2);
    REQUIRE(result.genetic().joint().probabilities.size() == 4);
    REQUIRE(
        result.genetic()
            .get<gelex::GeneticMode::D>()
            .annotation_coefficients.size()
        == 2);
    REQUIRE(result.total_genetic_variance().explained_variance.mean > 0.0);
}
