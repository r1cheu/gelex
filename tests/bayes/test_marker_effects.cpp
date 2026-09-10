// Copyright 2026 RuLei Chen
// SPDX-License-Identifier: Apache-2.0

#include <Eigen/Core>
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <cstddef>
#include <filesystem>
#include <fstream>
#include <string>
#include <type_traits>
#include <vector>

#include "gelex/bayes/builtin_method.h"
#include "gelex/bayes/draws.h"
#include "gelex/bayes/genetic/marker_effect_traits.h"
#include "gelex/bayes/genotype/design.h"
#include "gelex/bayes/marker_effects.h"
#include "gelex/bayes/prior.h"
#include "gelex/bayes/state.h"
#include "gelex/bayes/variance/budget.h"
#include "gelex/data/reader.h"
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

auto names_of(const gelex::MarkerEffectTable& table) -> std::vector<std::string>
{
    return {table.names().begin(), table.names().end()};
}

// var(x_A beta_A + x_D beta_D) of one marker, built column by column.
auto direct_marker_variance(
    const gelex::bayes::GeneticDesign& design,
    Eigen::Index marker,
    double beta_a,
    double beta_d) -> double
{
    Eigen::VectorXd total(design.rows());
    design.projection(gelex::GeneticMode::A).multiply(marker, beta_a, total);
    design.projection(gelex::GeneticMode::D).axpy(marker, beta_d, total);
    return gelex::vecvar(total, gelex::VarNormType::Population);
}

}  // namespace

TEST_CASE(
    "Marker effects scale posterior means by the genetic plus residual "
    "variance",
    "[bayes][marker_effects]")
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
            const auto value = static_cast<double>(draw + 1);
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

    const auto table = gelex::read_marker_effects<prior_type>(path, model);
    const gelex::DenseReader reader{path};
    const Eigen::VectorXd beta_a
        = reader.to_map<double>("genetic/A/coefficients").rowwise().mean();
    const Eigen::VectorXd beta_d
        = reader.to_map<double>("genetic/D/coefficients").rowwise().mean();
    const double residual
        = reader.to_map<double>("residual/variance").row(0).mean();
    const auto& design = model.genetic();

    Eigen::VectorXd gebv = Eigen::VectorXd::Zero(design.rows());
    for (Eigen::Index marker = 0; marker < design.cols(); ++marker)
    {
        design.projection(gelex::GeneticMode::A)
            .axpy(marker, beta_a(marker), gebv);
        design.projection(gelex::GeneticMode::D)
            .axpy(marker, beta_d(marker), gebv);
    }
    const double denominator
        = gelex::vecvar(gebv, gelex::VarNormType::Population) + residual;

    SECTION("columns follow the mode order and end with the total PVE")
    {
        REQUIRE(
            names_of(table)
            == std::vector<std::string>{
                "BETA_A", "SE_A", "PVE_A", "BETA_D", "SE_D", "PVE_D", "PVE"});
        REQUIRE(table.rows() == design.cols());
    }

    SECTION("BETA and PVE per mode use posterior means")
    {
        REQUIRE(table.column("BETA_A").isApprox(beta_a));
        REQUIRE(table.column("BETA_D").isApprox(beta_d));
        const Eigen::VectorXd expected_pve_a
            = beta_a.array().square()
              * design.projection(gelex::GeneticMode::A)
                    .col_var()
                    .transpose()
                    .array()
              / denominator;
        REQUIRE(table.column("PVE_A").isApprox(expected_pve_a));
        REQUIRE(table.column("PVE_A")(0) > 0.0);
    }

    SECTION("total PVE equals the direct variance of the A + D contribution")
    {
        const auto& pve = table.column("PVE");
        for (Eigen::Index marker = 0; marker < design.cols(); ++marker)
        {
            REQUIRE_THAT(
                pve(marker),
                WithinAbs(
                    direct_marker_variance(
                        design, marker, beta_a(marker), beta_d(marker))
                        / denominator,
                    1e-12));
        }
        REQUIRE(pve(1) != table.column("PVE_A")(1) + table.column("PVE_D")(1));
    }

    SECTION("write_marker_effects emits metadata followed by every column")
    {
        const auto out = (fixture.get_test_dir() / "rr.snpeff").string();
        gelex::write_marker_effects(out, design, table);

        std::ifstream file{out};
        std::vector<std::string> lines;
        for (std::string line; std::getline(file, line);)
        {
            lines.push_back(line);
        }
        REQUIRE(lines.size() == static_cast<std::size_t>(design.cols()) + 1);
        REQUIRE(
            lines[0]
            == "CHR\tSNP\tBP\tA1\tA2\tA1FREQ\tBETA_A\tSE_A\tPVE_A\tBETA_D\t"
               "SE_D\tPVE_D\tPVE");
        const auto keys = design.marker_metadata().index().keys();
        REQUIRE(lines[1].contains("\t" + keys[0] + "\t"));

        const auto read_back = gelex::read_snp_effects(out);
        REQUIRE(read_back.rows() == static_cast<std::size_t>(design.cols()));
        REQUIRE(read_back.contains("BETA_A"));
        REQUIRE(read_back.contains("BETA_D"));
        REQUIRE(read_back["BETA_A"].to_map<double>().isApprox(beta_a, 1e-8));

        gelex::MarkerEffectTable short_table{design.cols() - 1};
        REQUIRE_THROWS_AS(
            gelex::write_marker_effects(out, design, short_table),
            gelex::GelexException);
    }
}

TEST_CASE(
    "Marker effects of the joint family report PIP per mode and overall",
    "[bayes][marker_effects]")
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
            const auto value = static_cast<double>(draw + 1);
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

    const auto table = gelex::read_marker_effects<prior_type>(path, model);
    REQUIRE(
        names_of(table)
        == std::vector<std::string>{
            "BETA_A",
            "SE_A",
            "PVE_A",
            "PIP_A",
            "BETA_D",
            "SE_D",
            "PVE_D",
            "PIP_D",
            "PIP",
            "PVE"});
    REQUIRE(table.column("PIP_A").isApprox(Eigen::VectorXd{{1.0, 0.0}}));
    REQUIRE(table.column("PIP_D").isApprox(Eigen::VectorXd{{0.0, 1.0}}));
    REQUIRE(table.column("PIP").isApprox(Eigen::VectorXd{{1.0, 1.0}}));
    REQUIRE(table.column("BETA_A")(0) == 2.5);
}
