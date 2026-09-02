/*
 * Copyright 2026 RuLei Chen
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 * http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

#include <Eigen/Core>
#include <array>
#include <catch2/catch_test_macros.hpp>
#include <filesystem>
#include <fstream>
#include <iterator>
#include <string>
#include <utility>
#include <vector>

#include "gelex/bayes/draws.h"
#include "gelex/bayes/genetic_family.h"
#include "gelex/bayes/model.h"
#include "gelex/bayes/prior.h"
#include "gelex/bayes/recipe.h"
#include "gelex/bayes/result.h"
#include "gelex/bayes/result_io.h"
#include "gelex/bayes/state.h"
#include "gelex/bayes/variance_budget.h"
#include "gelex/data/covariates.h"
#include "gelex/data/dataframe/constants.h"
#include "gelex/data/fixed_design.h"
#include "gelex/genetic_mode.h"

#include "compact_genotype_fixture.h"
#include "file_fixture.h"
#include "random_design_fixture.h"

namespace
{

constexpr auto mode_a = gelex::GeneticModeSet{gelex::GeneticMode::A};
constexpr auto mode_ad = gelex::GeneticMode::A | gelex::GeneticMode::D;

auto read_text(const std::filesystem::path& path) -> std::string
{
    std::ifstream input{path};
    return std::string{
        std::istreambuf_iterator<char>{input},
        std::istreambuf_iterator<char>{}};
}

auto summary_row_keys(const std::filesystem::path& path)
    -> std::vector<std::string>
{
    std::ifstream input{path};
    std::string line;
    std::getline(input, line);
    REQUIRE(line == "term\tindex\tmean\tstddev");

    std::vector<std::string> keys;
    while (std::getline(input, line))
    {
        const auto first_separator = line.find('\t');
        const auto second_separator = line.find('\t', first_separator + 1);
        REQUIRE(first_separator != std::string::npos);
        REQUIRE(second_separator != std::string::npos);
        keys.push_back(line.substr(0, second_separator));
    }
    return keys;
}

auto make_parameter_model() -> gelex::BayesModel
{
    const auto discrete_column = std::string{"group"} + gelex::separator + "B";
    auto fixed = gelex::FixedDesign::make(
        gelex::QuantitativeCovariate{
            .names = {"shared"}, .X = Eigen::MatrixXd{{1.0}, {2.0}, {3.0}}},
        gelex::DiscreteCovariate{
            .terms
            = {{.name = "group", .levels = {"A", "B"}, .reference_level = "A"}},
            .X = Eigen::MatrixXd{{0.0}, {1.0}, {0.0}}});

    const std::array random_column_names{
        std::string{"shared"}, discrete_column};
    std::vector<gelex::bayes::RandomDesign> random;
    random.push_back(
        gelex::test::make_random_design(
            "batch",
            random_column_names,
            Eigen::MatrixXd{{1.0, 0.0}, {0.0, 1.0}, {1.0, 0.0}}));

    return gelex::BayesModel{
        Eigen::VectorXd{{1.0, 2.0, 3.0}},
        std::move(fixed),
        std::move(random),
        gelex::test::make_genetic_design(
            Eigen::MatrixXd{{0.0, 1.0}, {1.0, 0.0}, {2.0, 1.0}}, mode_a)};
}

auto collect_parameter_result(const std::filesystem::path& path)
{
    using Family = gelex::GaussianFamily<gelex::VarianceLayout::Pooled>;

    auto model = make_parameter_model();
    const auto prior = gelex::make_prior(
        gelex::BayesRecipe<mode_a, Family>{
            gelex::VarianceBudget{{.additive = 0.4, .random = 0.2}}},
        model);
    auto state = gelex::make_state(prior, model);
    state.fixed().coefficients = Eigen::VectorXd{{1.0, 2.0, 3.0}};
    state.random()[0].coefficients = Eigen::VectorXd{{4.0, 5.0}};
    state.random()[0].variance = 6.0;
    state.random()[0].fitted_values = Eigen::VectorXd{{0.0, 3.0, 0.0}};
    state.genetic().get<gelex::GeneticMode::A>().family_state.variance = 0.5;
    state.residual().variance = 4.0;

    gelex::BayesDraws draws(prior, model, path.string(), 2);
    draws.append(state);
    draws.append(state);
    return gelex::make_result(model, draws);
}

template <gelex::GeneticModeSet Modes, typename Family, typename Configure>
auto collect_genetic_result(
    const std::filesystem::path& path,
    Configure configure)
{
    auto model = gelex::test::make_compact_model(
        Eigen::MatrixXd{{0.0, 1.0}, {1.0, 0.0}, {2.0, 1.0}},
        Eigen::VectorXd{{1.0, 2.0, 3.0}},
        Modes);
    const auto prior = gelex::make_prior(
        gelex::BayesRecipe<Modes, Family>::defaults(), model);
    auto state = gelex::make_state(prior, model);
    configure(state);

    gelex::BayesDraws draws(prior, model, path.string(), 2);
    draws.append(state);
    draws.append(state);
    return gelex::make_result(model, draws);
}

}  // namespace

TEST_CASE(
    "Bayes parameter writer preserves result identifiers and coefficient "
    "order",
    "[bayes][result_io]")
{
    gelex::test::FileFixture fixture;
    const auto result = collect_parameter_result(
        fixture.get_test_dir() / "parameter_writer.draws");
    const auto prefix = (fixture.get_test_dir() / "opaque").string();

    gelex::write_params(result, prefix);

    const auto discrete_column = std::string{"group"} + gelex::separator + "B";
    const auto expected
        = std::string{
              "term\tindex\tcolumn_name\tmean\tstddev\n"
              "fixed/coefficients\t0\tIntercept\t1.00000000e+00\t"
              "0.00000000e+00\n"
              "fixed/coefficients\t1\tshared\t2.00000000e+00\t"
              "0.00000000e+00\n"
              "fixed/coefficients\t2\t"}
        + discrete_column
        + "\t3.00000000e+00\t0.00000000e+00\n"
          "random/batch/coefficients\t0\tshared\t4.00000000e+00\t"
          "0.00000000e+00\n"
          "random/batch/coefficients\t1\t"
        + discrete_column + "\t5.00000000e+00\t0.00000000e+00\n";
    REQUIRE(read_text(prefix + ".params") == expected);
    REQUIRE_FALSE(std::filesystem::exists(prefix));
}

TEST_CASE(
    "Bayes summary writer follows binary reservation order and separates "
    "parameters",
    "[bayes][result_io]")
{
    gelex::test::FileFixture fixture;
    const auto result = collect_parameter_result(
        fixture.get_test_dir() / "summary_writer.draws");
    const auto prefix = (fixture.get_test_dir() / "summary").string();

    gelex::write_summary(result, prefix);

    const auto expected = std::string{
        "term\tindex\tmean\tstddev\n"
        "random/batch/variance\t0\t6.00000000e+00\t0.00000000e+00\n"
        "random/batch/explained_variance\t0\t2.00000000e+00\t"
        "0.00000000e+00\n"
        "genetic/A/variance\t0\t5.00000000e-01\t0.00000000e+00\n"
        "residual/variance\t0\t4.00000000e+00\t0.00000000e+00\n"
        "genetic/A/explained_variance\t0\t0.00000000e+00\t0.00000000e+00\n"
        "genetic/A/heritability\t0\t0.00000000e+00\t0.00000000e+00\n"
        "genetic/total/explained_variance\t0\t0.00000000e+00\t"
        "0.00000000e+00\n"
        "genetic/total/heritability\t0\t0.00000000e+00\t0.00000000e+00\n"};
    const auto content = read_text(prefix + ".summary");
    REQUIRE(content == expected);
    REQUIRE_FALSE(content.contains("coefficients"));
    REQUIRE_FALSE(std::filesystem::exists(prefix));
}

TEST_CASE(
    "Bayes summary writer expands vectors and omits marker-sized fields",
    "[bayes][result_io]")
{
    gelex::test::FileFixture fixture;

    SECTION("vector posterior")
    {
        using Family = gelex::ScaledMixtureFamily<>;
        const auto result = collect_genetic_result<mode_a, Family>(
            fixture.get_test_dir() / "vector_summary.draws",
            [](auto& state)
            {
                auto& family = state.genetic()
                                   .template get<gelex::GeneticMode::A>()
                                   .family_state;
                family.variance = 0.8;
                family.probabilities = {0.5, 0.2, 0.15, 0.1, 0.05};
                family.fitted_values = Eigen::MatrixXd{
                    {0.0, 1.0, 2.0, 0.0},
                    {3.0, 1.0, 0.0, 0.0},
                    {0.0, 1.0, 4.0, 0.0}};
            });
        const auto prefix = (fixture.get_test_dir() / "vector").string();

        gelex::write_summary(result, prefix);
        const auto output = std::filesystem::path{prefix + ".summary"};

        REQUIRE(
            summary_row_keys(output)
            == std::vector<std::string>{
                "genetic/A/variance\t0",
                "genetic/A/probabilities\t0",
                "genetic/A/probabilities\t1",
                "genetic/A/probabilities\t2",
                "genetic/A/probabilities\t3",
                "genetic/A/probabilities\t4",
                "genetic/A/component_explained_variance\t0",
                "genetic/A/component_explained_variance\t1",
                "genetic/A/component_explained_variance\t2",
                "genetic/A/component_explained_variance\t3",
                "residual/variance\t0",
                "genetic/A/explained_variance\t0",
                "genetic/A/heritability\t0",
                "genetic/total/explained_variance\t0",
                "genetic/total/heritability\t0"});
        const auto content = read_text(output);
        REQUIRE_FALSE(content.contains("coefficients"));
        REQUIRE_FALSE(content.contains("assignment"));
    }

    SECTION("unpooled marker variance")
    {
        using Family = gelex::SpikeSlabFamily<
            gelex::VarianceLayout::Unpooled,
            gelex::UpdatePolicy::Sampled>;
        const auto result = collect_genetic_result<mode_a, Family>(
            fixture.get_test_dir() / "unpooled_summary.draws",
            [](auto& state)
            {
                auto& family = state.genetic()
                                   .template get<gelex::GeneticMode::A>()
                                   .family_state;
                family.variance = Eigen::VectorXd{{0.1, 0.2}};
                family.probability = 0.35;
            });
        const auto prefix = (fixture.get_test_dir() / "unpooled").string();

        gelex::write_summary(result, prefix);
        const auto output = std::filesystem::path{prefix + ".summary"};

        REQUIRE(
            summary_row_keys(output)
            == std::vector<std::string>{
                "genetic/A/probability\t0",
                "residual/variance\t0",
                "genetic/A/explained_variance\t0",
                "genetic/A/heritability\t0",
                "genetic/total/explained_variance\t0",
                "genetic/total/heritability\t0"});
        REQUIRE_FALSE(read_text(output).contains("genetic/A/variance"));
    }

    SECTION("joint family")
    {
        using Family = gelex::JointSpikeSlabFamily<>;
        const auto result = collect_genetic_result<mode_ad, Family>(
            fixture.get_test_dir() / "joint_summary.draws",
            [](auto& state)
            {
                auto& dominance = state.genetic()
                                      .template get<gelex::GeneticMode::D>()
                                      .family_state;
                dominance.variance = 0.25;
                dominance.positive_probability = 0.7;
                auto& joint = state.genetic().joint();
                joint.probabilities = {0.7, 0.1, 0.1, 0.1};
                joint.fitted_values = Eigen::MatrixXd{
                    {1.0, 0.0, 2.0, 0.0},
                    {2.0, 0.0, 2.0, 4.0},
                    {3.0, 3.0, 2.0, 2.0}};
            });
        const auto prefix = (fixture.get_test_dir() / "joint").string();

        gelex::write_summary(result, prefix);
        const auto output = std::filesystem::path{prefix + ".summary"};

        REQUIRE(
            summary_row_keys(output)
            == std::vector<std::string>{
                "genetic/A/variance\t0",
                "genetic/D/variance\t0",
                "genetic/D/positive_probability\t0",
                "genetic/joint/probabilities\t0",
                "genetic/joint/probabilities\t1",
                "genetic/joint/probabilities\t2",
                "genetic/joint/probabilities\t3",
                "genetic/joint/component_explained_variance\t0",
                "genetic/joint/component_explained_variance\t1",
                "genetic/joint/component_explained_variance\t2",
                "genetic/joint/component_explained_variance\t3",
                "residual/variance\t0",
                "genetic/A/explained_variance\t0",
                "genetic/D/explained_variance\t0",
                "genetic/A/heritability\t0",
                "genetic/D/heritability\t0",
                "genetic/total/explained_variance\t0",
                "genetic/total/heritability\t0"});
    }
}
