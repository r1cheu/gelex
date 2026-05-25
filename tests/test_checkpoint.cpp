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

#include <filesystem>
#include <memory>
#include <random>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

#include <Eigen/Core>
#include <catch2/catch_test_macros.hpp>

#include "gelex/algo/infer/mcmc/recipes.h"
#include "gelex/algo/infer/mcmc/solver.h"
#include "gelex/algo/infer/params.h"
#include "gelex/data/genotype/genotype.h"
#include "gelex/io/binary_writer.h"
#include "gelex/io/mcmc/checkpoint_reader.h"
#include "gelex/io/mcmc/checkpoint_writer.h"
#include "gelex/model/bayes/gaussian_prior.h"
#include "gelex/model/bayes/model.h"
#include "gelex/model/bayes/prior.h"
#include "gelex/model/bayes/prior_parameters.h"
#include "gelex/model/bayes/state.h"
#include "gelex/model/bayes/state_capabilities.h"
#include "gelex/types/fixed_effects.h"
#include "gelex/types/genetic_effect_type.h"
#include "genotype_fixture.h"

using namespace gelex;            // NOLINT
using namespace gelex::genotype;  // NOLINT

namespace
{

constexpr Eigen::Index kNumIndividuals = 10;
constexpr Eigen::Index kNumMarkers = 5;

auto make_genotype(Eigen::Index n_samples, Eigen::Index n_snps) -> Genotype
{
    Eigen::MatrixXd data = Eigen::MatrixXd::Random(n_samples, n_snps);
    auto mean = data.colwise().mean().transpose().eval();
    auto stddev = Eigen::VectorXd::Ones(n_snps);
    return test::GenotypeBuilder::build(
        std::move(data), std::move(mean), std::move(stddev));
}

auto make_model(Eigen::Index n_samples, Eigen::Index n_snps) -> BayesModel
{
    auto phenotype = Eigen::VectorXd::Random(n_samples);
    auto fixed = FixedEffect::build(n_samples);

    std::vector<bayes::GeneticEffect> genetics;
    genetics.emplace_back(GeneticMode::A, make_genotype(n_samples, n_snps));

    return BayesModel{phenotype, std::move(fixed), {}, std::move(genetics)};
}

auto make_variance(double initial_value) -> bayes::VarianceParameter
{
    return bayes::VarianceParameter(
        initial_value, bayes::ScaledInvChiSqPrior{4.0, 1.0});
}

auto make_random_prior(double initial_value) -> bayes::RandomPrior
{
    return bayes::RandomPrior{make_variance(initial_value)};
}

auto make_residual_prior(double initial_value) -> bayes::ResidualPrior
{
    return bayes::ResidualPrior{make_variance(initial_value)};
}

auto make_bayes_a_prior(Eigen::Index /*n_snps*/ = kNumMarkers)
    -> bayes::BayesPrior
{
    std::vector<std::unique_ptr<bayes::GeneticPrior>> genetics;
    genetics.push_back(
        std::make_unique<bayes::GaussianPrior>(
            GeneticMode::A,
            bayes::MarkerVariance{
                bayes::MarkerVarianceLayout::per_marker, make_variance(0.2)}));
    return bayes::BayesPrior(
        make_random_prior(0.5), std::move(genetics), make_residual_prior(1.0));
}

auto make_spike_slab_prior() -> bayes::BayesPrior
{
    std::vector<std::unique_ptr<bayes::GeneticPrior>> genetics;
    genetics.push_back(
        std::make_unique<bayes::SpikeSlabGaussianPrior>(
            GeneticMode::A,
            bayes::MarkerVariance{
                bayes::MarkerVarianceLayout::per_marker, make_variance(0.2)},
            bayes::MixtureProportion{
                bayes::SimplexParameter{
                    Eigen::VectorXd{{0.8, 0.2}},
                    bayes::DirichletPrior{Eigen::VectorXd{{1.0, 1.0}}}},
                bayes::UpdatePolicy::sampled}));
    return bayes::BayesPrior(
        make_random_prior(0.5), std::move(genetics), make_residual_prior(1.0));
}

auto run_bayes_a(
    const BayesModel& model,
    Eigen::Index n_iters,
    std::string_view prefix,
    Eigen::Index seed) -> std::string
{
    mcmc::Params params{n_iters, 0, 1, n_iters};
    mcmc::Solver solver(
        params,
        mcmc::make_bayes_a_chain<mcmc::GeneticShape::a_only>,
        std::string(prefix),
        std::string(prefix));
    solver.run(model, make_bayes_a_prior(), seed);
    return std::string(prefix) + ".ckpt";
}

auto resume_bayes_a(
    const BayesModel& model,
    Eigen::Index n_iters,
    std::string_view prefix,
    const std::string& checkpoint_path) -> std::string
{
    mcmc::Params params{n_iters, 0, 1, n_iters};
    mcmc::Solver solver(
        params,
        mcmc::make_bayes_a_chain<mcmc::GeneticShape::a_only>,
        std::string(prefix),
        std::string(prefix));
    solver.resume(model, make_bayes_a_prior(), checkpoint_path);
    return std::string(prefix) + ".ckpt";
}

auto read_bayes_state_checkpoint(
    const BayesModel& model,
    const std::string& path) -> BayesState
{
    auto prior = make_bayes_a_prior();
    BayesState state(model, prior);
    static_cast<void>(read_checkpoint(path, state));
    return state;
}

auto fill_state(BayesState& state) -> void
{
    state.fixed().coeffs
        = Eigen::VectorXd::LinSpaced(state.fixed().coeffs.size(), 1.0, 2.0);

    auto& genetic = *state.genetic(GeneticMode::A);
    genetic.coeffs
        = Eigen::VectorXd::LinSpaced(genetic.coeffs.size(), 0.1, 1.0);
    genetic.u = Eigen::VectorXd::LinSpaced(genetic.u.size(), -0.5, 0.5);
    genetic.variance = 0.3;
    genetic.heritability = 0.4;

    state.residual().variance = 0.7;
    state.residual().y_adj
        = Eigen::VectorXd::Constant(state.residual().y_adj.size(), 1.5);
}

}  // namespace

TEST_CASE("BayesState checkpoint resume matches continuous run", "[checkpoint]")
{
    const std::string prefix_cont = "/tmp/gelex_test_ckpt_cont";
    const std::string prefix_first = "/tmp/gelex_test_ckpt_first";
    const std::string prefix_resume = "/tmp/gelex_test_ckpt_resume";

    auto model = make_model(30, 8);

    auto ckpt_cont = run_bayes_a(model, 20, prefix_cont, 42);
    auto ckpt_first = run_bayes_a(model, 10, prefix_first, 42);
    auto ckpt_resume = resume_bayes_a(model, 10, prefix_resume, ckpt_first);

    auto state_cont = read_bayes_state_checkpoint(model, ckpt_cont);
    auto state_resume = read_bayes_state_checkpoint(model, ckpt_resume);

    CHECK(state_cont.fixed().coeffs.isApprox(state_resume.fixed().coeffs, 0.0));
    const auto& genetic_cont = *state_cont.genetic(GeneticMode::A);
    const auto& genetic_resume = *state_resume.genetic(GeneticMode::A);
    CHECK(genetic_cont.coeffs.isApprox(genetic_resume.coeffs, 0.0));
    CHECK(genetic_cont.u.isApprox(genetic_resume.u, 0.0));
    CHECK(genetic_cont.variance == genetic_resume.variance);
    CHECK(state_cont.residual().variance == state_resume.residual().variance);
    CHECK(state_cont.residual().y_adj.isApprox(
        state_resume.residual().y_adj, 0.0));

    std::filesystem::remove(ckpt_cont);
    std::filesystem::remove(ckpt_first);
    std::filesystem::remove(ckpt_resume);
}

TEST_CASE("BayesState checkpoint round-trip preserves records", "[checkpoint]")
{
    const std::string prefix = "/tmp/gelex_test_bayes_state_ckpt_roundtrip";
    const std::string ckpt_path = prefix + ".ckpt";

    auto model = make_model(kNumIndividuals, kNumMarkers);
    auto prior = make_spike_slab_prior();
    BayesState state(model, prior);
    fill_state(state);

    auto& proportion = state.genetic_block_for(GeneticMode::A)
                           ->prior_state()
                           .require<bayes::ProportionStateCap>()
                           .proportion()[0];
    proportion.assignment = Eigen::VectorXi::Ones(kNumMarkers);
    proportion.count = Eigen::VectorXi{{0, static_cast<int>(kNumMarkers)}};
    proportion.value = Eigen::VectorXd{{0.1, 0.9}};

    std::mt19937_64 rng(123);
    rng();
    write_checkpoint(state, rng, prefix);

    auto restored_prior = make_spike_slab_prior();
    BayesState restored(model, restored_prior);
    auto restored_rng = read_checkpoint(ckpt_path, restored);

    CHECK(restored.fixed().coeffs.isApprox(state.fixed().coeffs));
    const auto& genetic = *state.genetic(GeneticMode::A);
    const auto& restored_genetic = *restored.genetic(GeneticMode::A);
    CHECK(restored_genetic.coeffs.isApprox(genetic.coeffs));
    CHECK(restored_genetic.u.isApprox(genetic.u));
    CHECK(restored_genetic.variance == genetic.variance);
    CHECK(restored_genetic.heritability == genetic.heritability);
    const auto& restored_proportion
        = restored.genetic_block_for(GeneticMode::A)
              ->prior_state()
              .require<bayes::ProportionStateCap>()
              .proportion()[0];
    CHECK(restored_proportion.assignment.isApprox(proportion.assignment));
    CHECK(restored_proportion.count.isApprox(proportion.count));
    CHECK(restored_proportion.value.isApprox(proportion.value));
    CHECK(restored.residual().y_adj.isApprox(state.residual().y_adj));
    CHECK(restored.residual().variance == state.residual().variance);
    CHECK(restored_rng() == rng());

    std::filesystem::remove(ckpt_path);
}

TEST_CASE(
    "BayesState checkpoint atomic write leaves no tmp file",
    "[checkpoint]")
{
    const std::string prefix = "/tmp/gelex_test_ckpt_atomic";
    const std::string ckpt_path = prefix + ".ckpt";
    const std::string tmp_path = prefix + ".ckpt.tmp";

    auto model = make_model(kNumIndividuals, kNumMarkers);
    auto prior = make_bayes_a_prior();
    BayesState state(model, prior);
    fill_state(state);
    std::mt19937_64 rng(99);

    write_checkpoint(state, rng, prefix);

    REQUIRE(std::filesystem::exists(ckpt_path));
    REQUIRE_FALSE(std::filesystem::exists(tmp_path));

    std::filesystem::remove(ckpt_path);
}

TEST_CASE("legacy checkpoint is rejected by BayesState resume", "[checkpoint]")
{
    const std::string prefix = "/tmp/gelex_test_legacy_ckpt_reject";
    const std::string ckpt_path = prefix + ".ckpt";

    auto model = make_model(kNumIndividuals, kNumMarkers);
    {
        io::BinaryWriter writer(ckpt_path);
        const auto handle = writer.reserve<double>("legacy_marker", 1, 1);
        writer.write(handle, 1.0);
    }

    auto prior = make_bayes_a_prior();
    BayesState state(model, prior);
    REQUIRE_THROWS_AS(read_checkpoint(ckpt_path, state), GelexException);

    std::filesystem::remove(ckpt_path);
}

TEST_CASE(
    "BayesState resume throws on checkpoint dimension mismatch",
    "[checkpoint]")
{
    const std::string prefix = "/tmp/gelex_test_ckpt_dimcheck";
    const std::string ckpt_path = prefix + ".ckpt";

    auto model_write = make_model(kNumIndividuals, kNumMarkers);
    auto prior_write = make_bayes_a_prior();
    BayesState state(model_write, prior_write);
    fill_state(state);
    std::mt19937_64 rng(42);

    write_checkpoint(state, rng, prefix);

    auto model_mismatch = make_model(kNumIndividuals, kNumMarkers + 3);
    auto prior_read = make_bayes_a_prior();
    BayesState mismatch_state(model_mismatch, prior_read);
    REQUIRE_THROWS_AS(
        read_checkpoint(ckpt_path, mismatch_state), GelexException);

    std::filesystem::remove(ckpt_path);
}
