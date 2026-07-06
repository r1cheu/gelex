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
#include <array>
#include <cmath>
#include <string>
#include <string_view>
#include <utility>
#include <variant>
#include <vector>

#include <Eigen/Core>
#include <catch2/catch_test_macros.hpp>

#include "file_fixture.h"
#include "gelex/algo/mcmc/records.h"
#include "gelex/bayes/design.h"
#include "gelex/bayes/genetic/gaussian_prior.h"
#include "gelex/bayes/genetic/half_normal_prior.h"
#include "gelex/bayes/genetic/parameters.h"
#include "gelex/bayes/model.h"
#include "gelex/bayes/prior.h"
#include "gelex/bayes/state.h"
#include "gelex/data/genotype.h"
#include "gelex/io/binary_reader.h"
#include "gelex/types/fixed_designs.h"
#include "gelex/types/genetic_mode.h"
#include "genotype_fixture.h"

namespace
{

auto make_variance(double value) -> gelex::bayes::VarianceParameter
{
    return gelex::bayes::VarianceParameter{
        value, gelex::bayes::ScaledInvChiSqPrior{4.0, 1.0}};
}

auto make_genotype(Eigen::MatrixXd data) -> gelex::Genotype
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
        gelex::FixedDesign::make(3),
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
    const std::vector<gelex::RecordEntry>& entries,
    std::string_view path) -> const gelex::RecordResult&
{
    const auto it = std::ranges::find(
        entries, std::string{path}, &gelex::RecordEntry::path);
    REQUIRE(it != entries.end());
    return it->value;
}

}  // namespace

TEST_CASE("Records stores traced BayesState fields", "[mcmc][mcmc_records]")
{
    auto model = make_model();
    auto prior = make_prior();
    gelex::BayesState state(model, prior);
    gelex::Records records{2, ""};

    auto& block
        = std::get<gelex::bayes::SingleGeneticBlockState>(state.genetics()[0]);
    auto& prior_state
        = std::get<gelex::bayes::SingleScaledMixtureGaussianState>(
            block.prior_state());

    state.fixed().coeffs = Eigen::VectorXd{{1.0}};
    state.random()[0].variance = 2.0;
    state.residual().variance = 5.0;
    prior_state.proportion() = Eigen::VectorXd{{0.1, 0.2, 0.3, 0.4}};
    prior_state.assignment() = Eigen::VectorXi{{0, 1}};
    records.store(model, state);

    state.fixed().coeffs = Eigen::VectorXd{{3.0}};
    state.random()[0].variance = 4.0;
    state.residual().variance = 9.0;
    prior_state.assignment() = Eigen::VectorXi{{2, 3}};
    records.store(model, state);

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
        "state/genetic_0/single/A/genetic/coeffs",
        "state/genetic_0/single/A/genetic/variance",
        "state/genetic_0/single/A/genetic/heritability",
        "state/genetic_0/single/A/prior_state/"
        "scaled_mixture_gaussian/variance",
        "state/genetic_0/single/A/prior_state/"
        "scaled_mixture_gaussian/component/gebv_var",
        "state/genetic_0/single/A/prior_state/"
        "scaled_mixture_gaussian/mixture/proportion",
        "state/genetic_0/single/A/prior_state/"
        "scaled_mixture_gaussian/mixture/assignment",
        "state/residual/variance"};
    REQUIRE(paths == expected_paths);

    const auto fixed_entry = std::ranges::find(
        entries, std::string{"state/fixed/coeffs"}, &gelex::RecordEntry::path);
    REQUIRE(fixed_entry != entries.end());
    REQUIRE(fixed_entry->names);
    REQUIRE(*fixed_entry->names == std::vector<std::string>{"Intercept"});

    const auto random_coeffs_entry = std::ranges::find(
        entries,
        std::string{"state/random_0/coeffs"},
        &gelex::RecordEntry::path);
    REQUIRE(random_coeffs_entry != entries.end());
    REQUIRE(random_coeffs_entry->names);
    REQUIRE(
        *random_coeffs_entry->names
        == std::vector<std::string>{"batch_a", "batch_b"});

    const auto genetic_variance_entry = std::ranges::find(
        entries,
        std::string{"state/genetic_0/single/A/genetic/variance"},
        &gelex::RecordEntry::path);
    REQUIRE(genetic_variance_entry != entries.end());
    REQUIRE(genetic_variance_entry->names);
    REQUIRE(*genetic_variance_entry->names == std::vector<std::string>{"σ²"});

    const auto heritability_entry = std::ranges::find(
        entries,
        std::string{"state/genetic_0/single/A/genetic/heritability"},
        &gelex::RecordEntry::path);
    REQUIRE(heritability_entry != entries.end());
    REQUIRE(heritability_entry->names);
    REQUIRE(*heritability_entry->names == std::vector<std::string>{"h²"});

    const auto marker_variance_entry = std::ranges::find(
        entries,
        std::string{"state/genetic_0/single/A/prior_state/"
                    "scaled_mixture_gaussian/variance"},
        &gelex::RecordEntry::path);
    REQUIRE(marker_variance_entry != entries.end());
    REQUIRE(marker_variance_entry->names);
    REQUIRE(
        *marker_variance_entry->names == std::vector<std::string>{"σ²_marker"});

    const auto component_entry = std::ranges::find(
        entries,
        std::string{"state/genetic_0/single/A/prior_state/"
                    "scaled_mixture_gaussian/component/gebv_var"},
        &gelex::RecordEntry::path);
    REQUIRE(component_entry != entries.end());
    REQUIRE(component_entry->names);
    REQUIRE(
        *component_entry->names
        == std::vector<std::string>{
            "σ²_component[0]", "σ²_component[1]", "σ²_component[2]"});

    const auto proportion_entry = std::ranges::find(
        entries,
        std::string{"state/genetic_0/single/A/prior_state/"
                    "scaled_mixture_gaussian/mixture/proportion"},
        &gelex::RecordEntry::path);
    REQUIRE(proportion_entry != entries.end());
    REQUIRE(proportion_entry->names);
    REQUIRE(
        *proportion_entry->names
        == std::vector<std::string>{"π[0]", "π[1]", "π[2]", "π[3]"});

    const auto assignment_entry = std::ranges::find(
        entries,
        std::string{"state/genetic_0/single/A/prior_state/"
                    "scaled_mixture_gaussian/mixture/assignment"},
        &gelex::RecordEntry::path);
    REQUIRE(assignment_entry != entries.end());
    REQUIRE_FALSE(assignment_entry->names);

    const auto fixed = std::get<gelex::RunningStatsResult>(
        require_record(entries, "state/fixed/coeffs"));
    REQUIRE(fixed.mean.isApprox(Eigen::VectorXd{{2.0}}));
    REQUIRE(fixed.stddev.isApprox(Eigen::VectorXd{{std::sqrt(2.0)}}));

    const auto random = std::get<gelex::RunningStatsResult>(
        require_record(entries, "state/random_0/variance"));
    REQUIRE(random.mean.isApprox(Eigen::VectorXd{{3.0}}));
    REQUIRE(random.stddev.isApprox(Eigen::VectorXd{{std::sqrt(2.0)}}));

    const auto residual = std::get<gelex::RunningStatsResult>(
        require_record(entries, "state/residual/variance"));
    REQUIRE(residual.mean.isApprox(Eigen::VectorXd{{7.0}}));
    REQUIRE(residual.stddev.isApprox(Eigen::VectorXd{{std::sqrt(8.0)}}));

    const Eigen::MatrixXd expected_probabilities{
        {0.5, 0.0, 0.5, 0.0}, {0.0, 0.5, 0.0, 0.5}};
    const auto assignment = std::get<gelex::CategoryProbResult>(require_record(
        entries,
        "state/genetic_0/single/A/prior_state/"
        "scaled_mixture_gaussian/mixture/assignment"));
    REQUIRE(assignment.value.isApprox(expected_probabilities));
}

TEST_CASE(
    "Records stores traced dominance sign categories",
    "[mcmc][mcmc_records]")
{
    std::vector<gelex::bayes::GeneticDesign> genetics;
    genetics.emplace_back(
        gelex::GeneticMode::A,
        make_genotype(Eigen::MatrixXd{{0.0, 1.0}, {1.0, 0.0}, {2.0, 1.0}}));
    genetics.emplace_back(
        gelex::GeneticMode::D,
        make_genotype(Eigen::MatrixXd{{1.0, 0.0}, {0.0, 1.0}, {1.0, 2.0}}));
    gelex::BayesModel model{
        Eigen::VectorXd{{1.0, 2.0, 3.0}},
        gelex::FixedDesign::make(3),
        {},
        std::move(genetics)};

    std::vector<gelex::bayes::GeneticPrior> genetic_priors;
    genetic_priors.emplace_back(
        gelex::bayes::JointGeneticPrior{
            gelex::bayes::JointHalfNormalMixturePrior{
                gelex::bayes::JointSharedMarkerVariance{std::array{
                    gelex::bayes::SharedMarkerVariance{make_variance(0.1)},
                    gelex::bayes::SharedMarkerVariance{make_variance(0.2)}}},
                gelex::bayes::MixtureProportion{
                    Eigen::VectorXd{{0.7, 0.1, 0.1, 0.1}}},
                gelex::bayes::ProbabilityParameter{
                    0.6, gelex::bayes::BetaPrior{1.0, 1.0}}}});
    gelex::bayes::BayesPrior prior{
        gelex::bayes::RandomPrior{make_variance(0.3)},
        std::move(genetic_priors),
        gelex::bayes::ResidualPrior{make_variance(0.4)}};
    gelex::BayesState state(model, prior);
    gelex::Records records{2, ""};

    auto& block
        = std::get<gelex::bayes::JointGeneticBlockState>(state.genetics()[0]);
    auto& prior_state = std::get<gelex::bayes::JointHalfNormalMixtureState>(
        block.prior_state());

    prior_state.dominance_sign().sign = Eigen::VectorXi{{0, 1}};
    records.store(model, state);

    prior_state.dominance_sign().sign = Eigen::VectorXi{{1, 1}};
    records.store(model, state);

    auto entries = std::move(records).take_results();
    const Eigen::MatrixXd expected_probabilities{{0.5, 0.5}, {0.0, 1.0}};
    const auto sign = std::get<gelex::CategoryProbResult>(require_record(
        entries,
        "state/genetic_0/joint/prior_state/"
        "joint_half_normal_mixture/dominance_sign/sign"));
    REQUIRE(sign.value.isApprox(expected_probabilities));
}

TEST_CASE("Records handoff consumes stored results", "[mcmc][mcmc_records]")
{
    auto model = make_model();
    auto prior = make_prior();
    gelex::BayesState state(model, prior);
    gelex::Records records{1, ""};

    records.store(model, state);

    REQUIRE_FALSE(std::move(records).take_results().empty());
    REQUIRE(std::move(records).take_results().empty());
}

TEST_CASE("Records writes retained draws", "[mcmc][mcmc_records]")
{
    auto model = make_model();
    auto prior = make_prior();
    gelex::BayesState state(model, prior);

    gelex::test::FileFixture files;
    const auto draws_path = files.generate_random_file_path(".draws");
    gelex::Records records{2, draws_path.string()};

    auto& block
        = std::get<gelex::bayes::SingleGeneticBlockState>(state.genetics()[0]);
    auto& prior_state
        = std::get<gelex::bayes::SingleScaledMixtureGaussianState>(
            block.prior_state());

    state.fixed().coeffs = Eigen::VectorXd{{1.0}};
    state.random()[0].variance = 2.0;
    state.residual().variance = 5.0;
    prior_state.proportion() = Eigen::VectorXd{{0.1, 0.2, 0.3, 0.4}};
    prior_state.assignment() = Eigen::VectorXi{{0, 1}};
    records.store(model, state);

    state.fixed().coeffs = Eigen::VectorXd{{3.0}};
    state.random()[0].variance = 4.0;
    state.residual().variance = 9.0;
    prior_state.proportion() = Eigen::VectorXd{{0.4, 0.3, 0.2, 0.1}};
    prior_state.assignment() = Eigen::VectorXi{{2, 3}};
    records.store(model, state);

    const auto entries = std::move(records).take_results();
    gelex::BinaryReader reader(draws_path.string());
    REQUIRE(reader.n_sections() == entries.size());

    const auto fixed = reader.to_map<double>("state/fixed/coeffs");
    REQUIRE(fixed.rows() == 1);
    REQUIRE(fixed.cols() == 2);
    REQUIRE(fixed.isApprox(Eigen::MatrixXd{{1.0, 3.0}}));

    const auto random_variance
        = reader.to_map<double>("state/random_0/variance");
    REQUIRE(random_variance.isApprox(Eigen::MatrixXd{{2.0, 4.0}}));

    const auto residual = reader.to_map<double>("state/residual/variance");
    REQUIRE(residual.isApprox(Eigen::MatrixXd{{5.0, 9.0}}));

    const auto proportion = reader.to_map<double>(
        "state/genetic_0/single/A/prior_state/"
        "scaled_mixture_gaussian/mixture/proportion");
    REQUIRE(proportion.isApprox(
        Eigen::MatrixXd{{0.1, 0.4}, {0.2, 0.3}, {0.3, 0.2}, {0.4, 0.1}}));

    const auto assignment = reader.to_map<int>(
        "state/genetic_0/single/A/prior_state/"
        "scaled_mixture_gaussian/mixture/assignment");
    Eigen::Matrix<int, Eigen::Dynamic, Eigen::Dynamic> expected_assignment(
        2, 2);
    expected_assignment = Eigen::Matrix<int, 2, 2>{{0, 2}, {1, 3}};
    REQUIRE(assignment == expected_assignment);

    REQUIRE_FALSE(reader.contains("state/fixed/coeffs_names"));
    REQUIRE_FALSE(
        reader.contains("state/genetic_0/single/A/genetic/variance_name"));
}
