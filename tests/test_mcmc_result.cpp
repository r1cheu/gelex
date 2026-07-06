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

#include <array>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iterator>
#include <ranges>
#include <string>
#include <string_view>
#include <utility>
#include <variant>
#include <vector>

#include <Eigen/Core>
#include <catch2/catch_test_macros.hpp>

#include "file_fixture.h"
#include "gelex/algo/mcmc/records.h"
#include "gelex/algo/mcmc/result.h"
#include "gelex/bayes/design.h"
#include "gelex/bayes/genetic/gaussian_prior.h"
#include "gelex/bayes/genetic/parameters.h"
#include "gelex/bayes/model.h"
#include "gelex/bayes/prior.h"
#include "gelex/bayes/state.h"
#include "gelex/data/genotype.h"
#include "gelex/exception.h"
#include "gelex/io/mcmc.h"
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

auto make_records(const gelex::BayesModel& model, gelex::BayesState& state)
    -> gelex::Records
{
    gelex::Records records{2, ""};
    auto& block
        = std::get<gelex::bayes::SingleGeneticBlockState>(state.genetics()[0]);
    auto& prior_state
        = std::get<gelex::bayes::SingleScaledMixtureGaussianState>(
            block.prior_state());

    state.fixed().coeffs = Eigen::VectorXd{{1.0}};
    state.random()[0].coeffs = Eigen::VectorXd{{10.0, 20.0}};
    state.random()[0].variance = 2.0;
    block.state().coeffs = Eigen::VectorXd{{0.5, 1.0}};
    state.residual().variance = 5.0;
    prior_state.assignment() = Eigen::VectorXi{{0, 1}};
    records.store(model, state);

    state.fixed().coeffs = Eigen::VectorXd{{3.0}};
    state.random()[0].coeffs = Eigen::VectorXd{{14.0, 28.0}};
    state.random()[0].variance = 4.0;
    block.state().coeffs = Eigen::VectorXd{{1.5, 2.0}};
    state.residual().variance = 9.0;
    prior_state.assignment() = Eigen::VectorXi{{2, 3}};
    records.store(model, state);
    return records;
}

}  // namespace

TEST_CASE("Result owns finalized record values", "[mcmc][mcmc_result]")
{
    auto model = make_model();
    auto prior = make_prior();
    gelex::BayesState state(model, prior);
    auto records = make_records(model, state);

    gelex::Result result{std::move(records), model, 2};

    REQUIRE(result.samples_collected() == 2);
    REQUIRE(result.phenotype_variance() == model.phenotype_variance());
    REQUIRE_FALSE(result.records().empty());
    REQUIRE(result.records()[0].path == "state/fixed/coeffs");
    REQUIRE(result.records()[0].names);
    REQUIRE(
        *result.records()[0].names == std::vector<std::string>{"Intercept"});

    const auto& fixed
        = std::get<gelex::RunningStatsResult>(result.get("state/fixed/coeffs"));
    REQUIRE(fixed.mean.isApprox(Eigen::VectorXd{{2.0}}));
    REQUIRE(fixed.stddev.isApprox(Eigen::VectorXd{{std::sqrt(2.0)}}));

    const auto& residual = std::get<gelex::RunningStatsResult>(
        result.get("state/residual/variance"));
    REQUIRE(residual.mean.isApprox(Eigen::VectorXd{{7.0}}));
    REQUIRE(residual.stddev.isApprox(Eigen::VectorXd{{std::sqrt(8.0)}}));

    const auto& pve = std::get<gelex::RunningStatsResult>(
        result.get("state/genetic_0/single/A/genetic/pve"));
    REQUIRE(pve.mean.isApprox(Eigen::VectorXd{{1.0, 0.75}}));
    REQUIRE(pve.stddev.isApprox(Eigen::VectorXd::Zero(2)));
    const auto pve_record = std::ranges::find(
        result.records(),
        std::string{"state/genetic_0/single/A/genetic/pve"},
        &gelex::RecordEntry::path);
    REQUIRE(pve_record != result.records().end());
    REQUIRE_FALSE(pve_record->names);
    REQUIRE_THROWS_AS(result.get("state/fixed/pve"), gelex::GelexException);
    REQUIRE_THROWS_AS(result.get("state/random_0/pve"), gelex::GelexException);

    REQUIRE_THROWS_AS(result.get("state/genetic/pve"), gelex::GelexException);

    constexpr std::string_view assignment_path{
        "state/genetic_0/single/A/prior_state/"
        "scaled_mixture_gaussian/mixture/assignment"};
    const Eigen::MatrixXd expected_probabilities{
        {0.5, 0.0, 0.5, 0.0}, {0.0, 0.5, 0.0, 0.5}};
    const auto& assignment
        = std::get<gelex::CategoryProbResult>(result.get(assignment_path));
    const auto& assignment_again
        = std::get<gelex::CategoryProbResult>(result.get(assignment_path));
    REQUIRE(assignment.value.isApprox(expected_probabilities));
    REQUIRE(assignment_again.value.isApprox(expected_probabilities));

    const auto& pip = std::get<gelex::RunningStatsResult>(result.get(
        "state/genetic_0/single/A/prior_state/"
        "scaled_mixture_gaussian/mixture/pip"));
    REQUIRE(pip.mean.isApprox(Eigen::VectorXd{{0.5, 1.0}}));
    REQUIRE(pip.stddev.isApprox(Eigen::VectorXd::Zero(2)));
    const auto pip_record = std::ranges::find(
        result.records(),
        std::string{"state/genetic_0/single/A/prior_state/"
                    "scaled_mixture_gaussian/mixture/pip"},
        &gelex::RecordEntry::path);
    REQUIRE(pip_record != result.records().end());
    REQUIRE_FALSE(pip_record->names);

    const gelex::RecordEntry* assignment_record = nullptr;
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
    auto records = make_records(model, state);
    gelex::Result result{std::move(records), model, 2};

    REQUIRE_THROWS_AS(result.get("state/missing/path"), gelex::GelexException);
}

TEST_CASE("Result derives joint genetic PIP by effect", "[mcmc][mcmc_result]")
{
    std::vector<gelex::bayes::GeneticDesign> genetics;
    genetics.emplace_back(
        gelex::GeneticMode::A,
        make_genotype(
            Eigen::MatrixXd{
                {-1.0, 1.0 / 3.0}, {0.0, -2.0 / 3.0}, {1.0, 1.0 / 3.0}}));
    genetics.emplace_back(
        gelex::GeneticMode::D,
        make_genotype(
            Eigen::MatrixXd{
                {-1.0, 1.0 / 3.0}, {0.0, -2.0 / 3.0}, {1.0, 1.0 / 3.0}}));
    gelex::BayesModel model{
        Eigen::VectorXd{{1.0, 2.0, 3.0}},
        gelex::FixedDesign::make(3),
        {},
        std::move(genetics)};

    std::vector<gelex::bayes::GeneticPrior> priors;
    priors.emplace_back(
        gelex::bayes::JointGeneticPrior{gelex::bayes::JointGaussianMixturePrior{
            gelex::bayes::JointSharedMarkerVariance{std::array{
                gelex::bayes::SharedMarkerVariance{make_variance(0.2)},
                gelex::bayes::SharedMarkerVariance{make_variance(0.3)}}},
            gelex::bayes::MixtureProportion{
                Eigen::VectorXd{{0.25, 0.25, 0.25, 0.25}}}}});
    gelex::bayes::BayesPrior prior{
        gelex::bayes::RandomPrior{make_variance(0.4)},
        std::move(priors),
        gelex::bayes::ResidualPrior{make_variance(0.5)}};

    gelex::BayesState state(model, prior);
    auto& block
        = std::get<gelex::bayes::JointGeneticBlockState>(state.genetics()[0]);
    auto& prior_state = std::get<gelex::bayes::JointGaussianMixtureState>(
        block.prior_state());
    gelex::Records records{2, ""};

    block.state(gelex::GeneticMode::A).coeffs = Eigen::VectorXd{{0.5, 1.0}};
    block.state(gelex::GeneticMode::D).coeffs = Eigen::VectorXd{{1.0, 0.5}};
    prior_state.assignment() = Eigen::VectorXi{{0, 1}};
    records.store(model, state);
    block.state(gelex::GeneticMode::A).coeffs = Eigen::VectorXd{{1.5, 2.0}};
    block.state(gelex::GeneticMode::D).coeffs = Eigen::VectorXd{{0.5, 1.5}};
    prior_state.assignment() = Eigen::VectorXi{{2, 3}};
    records.store(model, state);

    gelex::Result result{std::move(records), model, 2};
    constexpr std::string_view assignment_path{
        "state/genetic_0/joint/prior_state/"
        "joint_mixture_gaussian/mixture/assignment"};
    const Eigen::MatrixXd expected_probabilities{
        {0.5, 0.0, 0.5, 0.0}, {0.0, 0.5, 0.0, 0.5}};
    const auto& assignment
        = std::get<gelex::CategoryProbResult>(result.get(assignment_path));
    REQUIRE(assignment.value.isApprox(expected_probabilities));

    const auto& additive_pip = std::get<gelex::RunningStatsResult>(result.get(
        "state/genetic_0/joint/prior_state/"
        "joint_mixture_gaussian/mixture/A/pip"));
    const auto& dominance_pip = std::get<gelex::RunningStatsResult>(result.get(
        "state/genetic_0/joint/prior_state/"
        "joint_mixture_gaussian/mixture/D/pip"));
    REQUIRE(additive_pip.mean.isApprox(Eigen::VectorXd{{0.0, 1.0}}));
    REQUIRE(dominance_pip.mean.isApprox(Eigen::VectorXd{{0.5, 0.5}}));
    REQUIRE(additive_pip.stddev.isApprox(Eigen::VectorXd::Zero(2)));
    REQUIRE(dominance_pip.stddev.isApprox(Eigen::VectorXd::Zero(2)));
    const auto& total_pve
        = std::get<gelex::RunningStatsResult>(result.get("state/genetic/pve"));
    REQUIRE(
        total_pve.mean.isApprox(Eigen::VectorXd{{3.0625, 2.0833333333333335}}));
    REQUIRE_THROWS_AS(
        result.get(
            "state/genetic_0/joint/prior_state/"
            "joint_mixture_gaussian/mixture/pip"),
        gelex::GelexException);

    gelex::test::FileFixture files;
    auto bim_path = files.create_text_file(
        "1\trs1\t0\t100\tA\tG\n"
        "2\trs2\t0\t200\tC\tT\n",
        ".bim");
    auto prefix = files.get_test_dir() / "joint_mcmc_snp";
    gelex::write_snp_eff(result, model, bim_path, prefix.string());
    auto snp_path = prefix;
    snp_path += ".snpeff";
    REQUIRE(std::filesystem::exists(snp_path));
    std::ifstream input(snp_path);
    const std::string content{
        std::istreambuf_iterator<char>{input},
        std::istreambuf_iterator<char>{}};
    REQUIRE(
        content.find(
            "CHR\tSNP\tBP\tA1\tA2\tA1FREQ\tBETA_A\tSE_A\tPVE_A\t"
            "PIP_A\tBETA_D\tSE_D\tPVE_D\tPIP_D\tPVE\n")
        == 0);
}

TEST_CASE("write_summary writes user-facing summary", "[mcmc][mcmc_result]")
{
    auto model = make_model();
    auto prior = make_prior();
    gelex::BayesState state(model, prior);
    auto records = make_records(model, state);
    gelex::Result result{std::move(records), model, 2};

    gelex::test::FileFixture files;
    auto prefix = files.get_test_dir() / "mcmc_result";
    gelex::write_summary(result, prefix.string());

    auto summary_path = prefix;
    summary_path += ".summary";
    REQUIRE(std::filesystem::exists(summary_path));

    std::ifstream input(summary_path);
    const std::string content{
        std::istreambuf_iterator<char>{input},
        std::istreambuf_iterator<char>{}};

    REQUIRE(content.find("term\teffect\tmean\tstddev\n") == 0);
    REQUIRE(
        content.find("Intercept\t-\t2.00000000e+00\t") != std::string::npos);
    REQUIRE(content.find("σ²_e\t-\t7.00000000e+00\t") != std::string::npos);
    REQUIRE(content.find("σ²\tA\t") != std::string::npos);
    REQUIRE(content.find("σ²_marker\tA\t") != std::string::npos);
    REQUIRE(content.find("π[0]\tA\t") != std::string::npos);
    REQUIRE(content.find("additive") == std::string::npos);
    REQUIRE(content.find("dominance") == std::string::npos);
    REQUIRE(content.find("state/") == std::string::npos);
    REQUIRE(content.find("assignment") == std::string::npos);
    REQUIRE(content.find("coeffs") == std::string::npos);
    REQUIRE(content.find("pve") == std::string::npos);
    REQUIRE(content.find("pip") == std::string::npos);
}

TEST_CASE("write_params writes fixed and random effects", "[mcmc][mcmc_result]")
{
    auto model = make_model();
    auto prior = make_prior();
    gelex::BayesState state(model, prior);
    auto records = make_records(model, state);
    gelex::Result result{std::move(records), model, 2};

    gelex::test::FileFixture files;
    auto prefix = files.get_test_dir() / "mcmc_params";
    gelex::write_params(result, prefix.string());

    auto params_path = prefix;
    params_path += ".param";
    REQUIRE(std::filesystem::exists(params_path));

    std::ifstream input(params_path);
    const std::string content{
        std::istreambuf_iterator<char>{input},
        std::istreambuf_iterator<char>{}};

    REQUIRE(content.find("term\tmean\tstddev\n") == 0);
    REQUIRE(content.find("Intercept\t2.00000000e+00\t") != std::string::npos);
    REQUIRE(content.find("batch_a\t1.20000000e+01\t") != std::string::npos);
    REQUIRE(content.find("batch_b\t2.40000000e+01\t") != std::string::npos);
    REQUIRE(content.find("σ²") == std::string::npos);
    REQUIRE(content.find("π[") == std::string::npos);
    REQUIRE(content.find("state/") == std::string::npos);
    REQUIRE(content.find("assignment") == std::string::npos);
}

TEST_CASE(
    "write_snp_eff writes dynamic SNP effect columns",
    "[mcmc][mcmc_result]")
{
    auto model = make_model();
    auto prior = make_prior();
    gelex::BayesState state(model, prior);
    auto records = make_records(model, state);
    gelex::Result result{std::move(records), model, 2};

    gelex::test::FileFixture files;
    auto bim_path = files.create_text_file(
        "1\trs1\t0\t100\tA\tG\n"
        "2\trs2\t0\t200\tC\tT\n",
        ".bim");
    auto prefix = files.get_test_dir() / "mcmc_snp";
    gelex::write_snp_eff(result, model, bim_path, prefix.string());

    auto snp_path = prefix;
    snp_path += ".snpeff";
    REQUIRE(std::filesystem::exists(snp_path));

    std::ifstream input(snp_path);
    const std::string content{
        std::istreambuf_iterator<char>{input},
        std::istreambuf_iterator<char>{}};

    REQUIRE(
        content.find(
            "CHR\tSNP\tBP\tA1\tA2\tA1FREQ\tBETA_A\tSE_A\tPVE_A\t"
            "PIP_A\n")
        == 0);
    REQUIRE(
        content.find("1\trs1\t100\tA\tG\t5.00000000e-01\t1.00000000e+00\t")
        != std::string::npos);
    REQUIRE(
        content.find("\t1.00000000e+00\t5.00000000e-01\n")
        != std::string::npos);
    REQUIRE(
        content.find("2\trs2\t200\tC\tT\t3.33333333e-01\t1.50000000e+00\t")
        != std::string::npos);
    REQUIRE(content.find("BETA_D") == std::string::npos);
    REQUIRE(content.find("PIP_D") == std::string::npos);
    REQUIRE(content.find("PIP_A\tPVE") == std::string::npos);
}
