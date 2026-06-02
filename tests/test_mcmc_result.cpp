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

#include <cmath>
#include <filesystem>
#include <fstream>
#include <iterator>
#include <string>
#include <string_view>
#include <utility>
#include <variant>
#include <vector>

#include <Eigen/Core>
#include <catch2/catch_test_macros.hpp>

#include "gelex/algo/infer/mcmc/records.h"
#include "gelex/algo/infer/mcmc/result.h"
#include "gelex/bayes/design.h"
#include "gelex/bayes/genetic/gaussian_prior.h"
#include "gelex/bayes/genetic/parameters.h"
#include "gelex/bayes/model.h"
#include "gelex/bayes/prior.h"
#include "gelex/bayes/state.h"
#include "gelex/data/genotype/genotype.h"
#include "gelex/exception.h"
#include "gelex/io/mcmc/result_writer.h"
#include "gelex/types/fixed_designs.h"
#include "gelex/types/genetic_effect_type.h"
#include "file_fixture.h"
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
    std::vector<gelex::bayes::GeneticDesign> genetics;
    genetics.emplace_back(
        gelex::GeneticMode::A,
        make_genotype(Eigen::MatrixXd{{0.0, 1.0}, {1.0, 0.0}, {2.0, 1.0}}));

    return gelex::BayesModel{
        Eigen::VectorXd{{1.0, 2.0, 3.0}},
        gelex::FixedDesign::build(3),
        {},
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

auto make_records(
    const gelex::BayesModel& model,
    const gelex::bayes::BayesPrior& prior,
    gelex::BayesState& state)
    -> gelex::mcmc::Records
{
    gelex::mcmc::Records records;
    auto& block = std::get<gelex::bayes::SingleGeneticBlockState>(
        state.genetics()[0]);
    auto& prior_state
        = std::get<gelex::bayes::SingleScaledMixtureGaussianState>(
            block.prior_state());

    state.fixed().coeffs = Eigen::VectorXd{{1.0}};
    state.residual().variance = 5.0;
    prior_state.assignment() = Eigen::VectorXi{{0, 1}};
    records.store(model, prior, state);

    state.fixed().coeffs = Eigen::VectorXd{{3.0}};
    state.residual().variance = 9.0;
    prior_state.assignment() = Eigen::VectorXi{{2, 3}};
    records.store(model, prior, state);
    return records;
}

}  // namespace

TEST_CASE("Result owns finalized record values", "[mcmc][mcmc_result]")
{
    auto model = make_model();
    auto prior = make_prior();
    gelex::BayesState state(model, prior);
    auto records = make_records(model, prior, state);

    gelex::mcmc::Result result{std::move(records), model, 2};

    REQUIRE(result.samples_collected() == 2);
    REQUIRE(result.phenotype_variance() == model.phenotype_variance());
    REQUIRE_FALSE(result.records().empty());
    REQUIRE(result.records()[0].path == "state/fixed/coeffs");
    REQUIRE(result.records()[0].names);
    REQUIRE(*result.records()[0].names == std::vector<std::string>{"Intercept"});

    const auto& fixed = std::get<gelex::stats::RunningStatsResult>(
        result.get("state/fixed/coeffs"));
    REQUIRE(fixed.mean.isApprox(Eigen::VectorXd{{2.0}}));
    REQUIRE(fixed.stddev.isApprox(Eigen::VectorXd{{std::sqrt(2.0)}}));

    const auto& residual = std::get<gelex::stats::RunningStatsResult>(
        result.get("state/residual/variance"));
    REQUIRE(residual.mean.isApprox(Eigen::VectorXd{{7.0}}));
    REQUIRE(residual.stddev.isApprox(Eigen::VectorXd{{std::sqrt(8.0)}}));

    constexpr std::string_view assignment_path{
        "state/genetic_0/single/prior_state/"
        "scaled_mixture_gaussian/mixture/assignment"};
    const Eigen::MatrixXd expected_probabilities{
        {0.5, 0.0, 0.5, 0.0},
        {0.0, 0.5, 0.0, 0.5}};
    const auto& assignment = std::get<gelex::stats::CategoryProbResult>(
        result.get(assignment_path));
    const auto& assignment_again = std::get<gelex::stats::CategoryProbResult>(
        result.get(assignment_path));
    REQUIRE(assignment.value.isApprox(expected_probabilities));
    REQUIRE(assignment_again.value.isApprox(expected_probabilities));

    const gelex::mcmc::RecordEntry* assignment_record = nullptr;
    for (const auto& record : result.records())
    {
        if (record.path == assignment_path)
        {
            assignment_record = &record;
        }
    }
    REQUIRE(assignment_record != nullptr);
    REQUIRE_FALSE(assignment_record->names);
}

TEST_CASE("Result rejects missing record paths", "[mcmc][mcmc_result]")
{
    auto model = make_model();
    auto prior = make_prior();
    gelex::BayesState state(model, prior);
    auto records = make_records(model, prior, state);
    gelex::mcmc::Result result{std::move(records), model, 2};

    REQUIRE_THROWS_AS(result.get("state/missing/path"), gelex::GelexException);
}

TEST_CASE("write_result writes user-facing summary", "[mcmc][mcmc_result]")
{
    auto model = make_model();
    auto prior = make_prior();
    gelex::BayesState state(model, prior);
    auto records = make_records(model, prior, state);
    gelex::mcmc::Result result{std::move(records), model, 2};

    gelex::test::FileFixture files;
    auto prefix = files.get_test_dir() / "mcmc_result";
    gelex::mcmc::write_result(result, prefix);

    auto summary_path = prefix;
    summary_path += ".summary";
    REQUIRE(std::filesystem::exists(summary_path));

    std::ifstream input(summary_path);
    const std::string content{
        std::istreambuf_iterator<char>{input},
        std::istreambuf_iterator<char>{}};

    REQUIRE(content.find("term\tmean\tstddev\n") == 0);
    REQUIRE(content.find("Intercept\t2\t") != std::string::npos);
    REQUIRE(content.find("σ²_e\t7\t") != std::string::npos);
    REQUIRE(content.find("σ²_add\t") != std::string::npos);
    REQUIRE(content.find("π_add[0]\t") != std::string::npos);
    REQUIRE(content.find("state/") == std::string::npos);
    REQUIRE(content.find("assignment") == std::string::npos);
    REQUIRE(content.find("coeffs") == std::string::npos);
}
