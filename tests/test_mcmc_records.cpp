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

#include <algorithm>
#include <cmath>
#include <string>
#include <string_view>
#include <utility>
#include <variant>
#include <vector>

#include <Eigen/Core>
#include <catch2/catch_test_macros.hpp>

#include "gelex/algo/infer/mcmc/records.h"
#include "gelex/bayes/design.h"
#include "gelex/bayes/genetic/gaussian_prior.h"
#include "gelex/bayes/genetic/parameters.h"
#include "gelex/bayes/model.h"
#include "gelex/bayes/prior.h"
#include "gelex/bayes/state.h"
#include "gelex/data/genotype/genotype.h"
#include "gelex/types/fixed_designs.h"
#include "gelex/types/genetic_effect_type.h"
#include "genotype_fixture.h"

namespace
{

auto make_variance(double value) -> gelex::bayes::VarianceParameter
{
    return gelex::bayes::VarianceParameter{
        value, gelex::bayes::ScaledInvChiSqPrior{4.0, 1.0}};
}

auto make_genotype(Eigen::MatrixXd data) -> gelex::genotype::Genotype
{
    auto mean = data.colwise().mean().transpose().eval();
    auto stddev = Eigen::VectorXd::Ones(data.cols());
    return gelex::test::GenotypeBuilder::build(
        std::move(data), std::move(mean), std::move(stddev));
}

auto make_model() -> gelex::BayesModel
{
    std::vector<gelex::bayes::RandomDesign> random;
    random.emplace_back(
        "batch",
        std::vector<std::string>{"a", "b"},
        Eigen::MatrixXd{{1.0, 0.0}, {0.0, 1.0}, {1.0, 0.0}});

    std::vector<gelex::bayes::GeneticDesign> genetics;
    genetics.emplace_back(
        gelex::GeneticMode::A,
        make_genotype(Eigen::MatrixXd{{0.0, 1.0}, {1.0, 0.0}, {2.0, 1.0}}));

    return gelex::BayesModel{
        Eigen::VectorXd{{1.0, 2.0, 3.0}},
        gelex::FixedDesign::build(3),
        std::move(random),
        std::move(genetics)};
}

auto make_prior() -> gelex::bayes::BayesPrior
{
    std::vector<gelex::bayes::GeneticPrior> genetics;
    genetics.emplace_back(
        gelex::bayes::SingleGeneticPrior{
            gelex::bayes::SingleScaledMixtureGaussianPrior{
                gelex::GeneticMode::A,
                gelex::bayes::SharedMarkerVariance{make_variance(0.2)},
                Eigen::VectorXd{{0.0, 1.0, 2.0, 4.0}},
                gelex::bayes::MixtureProportion{
                    Eigen::VectorXd{{0.25, 0.25, 0.25, 0.25}}}}});

    return gelex::bayes::BayesPrior{
        gelex::bayes::RandomPrior{make_variance(0.3)},
        std::move(genetics),
        gelex::bayes::ResidualPrior{make_variance(0.4)}};
}

auto require_record(
    const std::vector<gelex::mcmc::RecordEntry>& entries,
    std::string_view path) -> const gelex::mcmc::RecordResult&
{
    const auto it = std::ranges::find(
        entries,
        std::string{path},
        &gelex::mcmc::RecordEntry::path);
    REQUIRE(it != entries.end());
    return it->value;
}

}  // namespace

TEST_CASE("Records stores traced BayesState fields", "[mcmc][mcmc_records]")
{
    auto model = make_model();
    auto prior = make_prior();
    gelex::BayesState state(model, prior);
    gelex::mcmc::Records records;

    auto& block = std::get<gelex::bayes::SingleGeneticBlockState>(
        state.genetics()[0]);
    auto& prior_state
        = std::get<gelex::bayes::SingleScaledMixtureGaussianState>(
            block.prior_state());

    state.fixed().coeffs = Eigen::VectorXd{{1.0}};
    state.random()[0].variance = 2.0;
    state.residual().variance = 5.0;
    prior_state.proportion() = Eigen::VectorXd{{0.1, 0.2, 0.3, 0.4}};
    prior_state.assignment() = Eigen::VectorXi{{0, 1}};
    records.store(state);

    state.fixed().coeffs = Eigen::VectorXd{{3.0}};
    state.random()[0].variance = 4.0;
    state.residual().variance = 9.0;
    prior_state.assignment() = Eigen::VectorXi{{2, 3}};
    records.store(state);

    auto entries = std::move(records).take_results();
    std::vector<std::string> paths;
    paths.reserve(entries.size());
    for (const auto& entry : entries)
    {
        paths.push_back(entry.path);
    }

    const std::vector<std::string> expected_paths{
        "state/fixed/coeffs",
        "state/random_0/coeffs",
        "state/random_0/variance",
        "state/genetic_0/single/genetic/coeffs",
        "state/genetic_0/single/genetic/variance",
        "state/genetic_0/single/genetic/heritability",
        "state/genetic_0/single/prior_state/"
        "scaled_mixture_gaussian/variance",
        "state/genetic_0/single/prior_state/"
        "scaled_mixture_gaussian/component/gebv_var",
        "state/genetic_0/single/prior_state/"
        "scaled_mixture_gaussian/mixture/proportion",
        "state/genetic_0/single/prior_state/"
        "scaled_mixture_gaussian/mixture/assignment",
        "state/residual/variance"};
    REQUIRE(paths == expected_paths);

    const auto fixed = std::get<gelex::stats::RunningStatsResult>(
        require_record(entries, "state/fixed/coeffs"));
    REQUIRE(fixed.mean.isApprox(Eigen::VectorXd{{2.0}}));
    REQUIRE(fixed.stddev.isApprox(Eigen::VectorXd{{std::sqrt(2.0)}}));

    const auto random = std::get<gelex::stats::RunningStatsResult>(
        require_record(entries, "state/random_0/variance"));
    REQUIRE(random.mean.isApprox(Eigen::VectorXd{{3.0}}));
    REQUIRE(random.stddev.isApprox(Eigen::VectorXd{{std::sqrt(2.0)}}));

    const auto residual = std::get<gelex::stats::RunningStatsResult>(
        require_record(entries, "state/residual/variance"));
    REQUIRE(residual.mean.isApprox(Eigen::VectorXd{{7.0}}));
    REQUIRE(residual.stddev.isApprox(Eigen::VectorXd{{std::sqrt(8.0)}}));

    const Eigen::MatrixXd expected_probabilities{
        {0.5, 0.0, 0.5, 0.0},
        {0.0, 0.5, 0.0, 0.5}};
    const auto assignment = std::get<gelex::stats::CategoryProbResult>(
        require_record(
            entries,
            "state/genetic_0/single/prior_state/"
            "scaled_mixture_gaussian/mixture/assignment"));
    REQUIRE(assignment.value.isApprox(expected_probabilities));
}

TEST_CASE("Records handoff consumes stored results", "[mcmc][mcmc_records]")
{
    auto model = make_model();
    auto prior = make_prior();
    gelex::BayesState state(model, prior);
    gelex::mcmc::Records records;

    records.store(state);

    REQUIRE_FALSE(std::move(records).take_results().empty());
    REQUIRE(std::move(records).take_results().empty());
}
