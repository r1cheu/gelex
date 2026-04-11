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
#include <random>
#include <string>
#include <utility>
#include <vector>

#include <Eigen/Core>
#include <catch2/catch_test_macros.hpp>

#include "gelex/algo/infer/mcmc.h"
#include "gelex/algo/infer/params.h"
#include "gelex/data/genotype/genotype_matrix.h"
#include "gelex/model/bayes/checkpoint.h"
#include "gelex/model/bayes/model.h"
#include "gelex/model/bayes/prior.h"
#include "gelex/model/bayes/prior_config.h"
#include "gelex/model/bayes/reader/checkpoint_reader.h"
#include "gelex/model/bayes/trait_model.h"
#include "gelex/model/bayes/writer/checkpoint_writer.h"
#include "gelex/types/bayes_method.h"
#include "gelex/types/fixed_effects.h"

using namespace gelex;  // NOLINT

namespace
{

constexpr Eigen::Index kNSamples = 10;
constexpr Eigen::Index kNSnps = 5;

auto make_geno_matrix(Eigen::Index n_samples, Eigen::Index n_snps)
    -> GenotypeMatrix
{
    auto data = Eigen::MatrixXd::Random(n_samples, n_snps);
    auto mean = data.colwise().mean().transpose().eval();
    auto stddev = Eigen::VectorXd::Ones(n_snps);
    return GenotypeMatrix(data, {}, std::move(mean), stddev);
}

auto make_bayes_a_model(Eigen::Index n_samples, Eigen::Index n_snps)
    -> std::pair<BayesModel, bayes::Priors>
{
    auto phenotype = Eigen::VectorXd::Random(n_samples);
    auto fixed = FixedEffect::build(n_samples);
    auto geno = make_geno_matrix(n_samples, n_snps);

    std::vector<bayes::GeneticEffect> genetics;
    genetics.emplace_back(
        GeneticMode::A, bayes::GenotypeStorage{std::move(geno)});

    double pheno_var
        = phenotype.array().square().mean() - std::pow(phenotype.mean(), 2.0);
    PriorSetConfig pc(BayesMethodConfig{BayesBase::A}, pheno_var);
    bayes::Priors priors(pc, genetics, 0);

    BayesModel model(phenotype, std::move(fixed), std::move(genetics));
    return {std::move(model), std::move(priors)};
}

auto run_bayes_a(
    BayesModel& model,
    const bayes::Priors& priors,
    Eigen::Index n_iters,
    std::string_view prefix,
    Eigen::Index seed) -> std::string
{
    MCMCParams params(n_iters, 0, 1, n_iters);
    MCMC mcmc(params, BayesA{}, std::string(prefix));
    mcmc.run(model, priors, seed, "");
    return std::string(prefix) + ".ckpt";
}

auto resume_bayes_a(
    BayesModel& model,
    Eigen::Index n_iters,
    std::string_view prefix,
    const std::string& ckpt_path) -> std::string
{
    auto ckpt = read_checkpoint(ckpt_path);
    MCMCParams params(n_iters, 0, 1, n_iters);
    MCMC mcmc(params, BayesA{}, std::string(prefix));
    mcmc.resume(model, std::move(ckpt), "");
    return std::string(prefix) + ".ckpt";
}

auto fill_state(BayesState& state) -> void
{
    state.fixed().coeffs
        = Eigen::VectorXd::LinSpaced(state.fixed().coeffs.size(), 1.0, 2.0);

    auto& gs = state.genetics()[0];
    gs.coeffs = Eigen::VectorXd::LinSpaced(gs.coeffs.size(), 0.1, 1.0);
    gs.u = Eigen::VectorXd::LinSpaced(gs.u.size(), -0.5, 0.5);
    gs.variance = 0.3;
    gs.heritability = 0.4;

    state.residual().variance = 0.7;
    state.residual().y_adj
        = Eigen::VectorXd::Constant(state.residual().y_adj.size(), 1.5);
}

}  // namespace

TEST_CASE("checkpoint resume produces bit-exact final state", "[checkpoint]")
{
    const std::string prefix_cont = "/tmp/gelex_test_ckpt_cont";
    const std::string prefix_first = "/tmp/gelex_test_ckpt_first";
    const std::string prefix_resume = "/tmp/gelex_test_ckpt_resume";

    Eigen::Index n_samples = 30;
    Eigen::Index n_snps = 8;
    auto [model, priors] = make_bayes_a_model(n_samples, n_snps);

    auto ckpt_cont = run_bayes_a(model, priors, 20, prefix_cont, 42);

    auto ckpt_first = run_bayes_a(model, priors, 10, prefix_first, 42);
    auto ckpt_resume = resume_bayes_a(model, 10, prefix_resume, ckpt_first);

    auto state_cont = read_checkpoint(ckpt_cont);
    auto state_resume = read_checkpoint(ckpt_resume);

    CHECK(state_cont.state.fixed().coeffs.isApprox(
        state_resume.state.fixed().coeffs, 0.0));

    const auto& gs_cont = state_cont.state.genetics()[0];
    const auto& gs_resume = state_resume.state.genetics()[0];
    CHECK(gs_cont.coeffs.isApprox(gs_resume.coeffs, 0.0));
    CHECK(gs_cont.u.isApprox(gs_resume.u, 0.0));
    CHECK(gs_cont.variance == gs_resume.variance);
    CHECK(gs_cont.marker_variance.isApprox(gs_resume.marker_variance, 0.0));

    // residual state
    CHECK(
        state_cont.state.residual().variance
        == state_resume.state.residual().variance);
    CHECK(state_cont.state.residual().y_adj.isApprox(
        state_resume.state.residual().y_adj, 0.0));

    CHECK(state_cont.rng() == state_resume.rng());

    std::filesystem::remove(ckpt_cont);
    std::filesystem::remove(ckpt_first);
    std::filesystem::remove(ckpt_resume);
}

TEST_CASE("checkpoint round-trip preserves all fields", "[checkpoint]")
{
    const std::string prefix = "/tmp/gelex_test_ckpt_roundtrip";
    const std::string ckpt_path = prefix + ".ckpt";

    auto [model, priors] = make_bayes_a_model(kNSamples, kNSnps);
    BayesState state(model, priors);
    fill_state(state);

    std::mt19937_64 rng(12345);
    rng();
    rng();
    rng();

    write_checkpoint(state, rng, priors, prefix);

    auto ckpt = read_checkpoint(ckpt_path);

    // fixed coeffs
    CHECK(ckpt.state.fixed().coeffs.isApprox(state.fixed().coeffs));

    // genetic state
    const auto& gs_orig = state.genetics()[0];
    const auto& gs_read = ckpt.state.genetics()[0];
    CHECK(gs_read.coeffs.isApprox(gs_orig.coeffs));
    CHECK(gs_read.u.isApprox(gs_orig.u));
    CHECK(gs_read.variance == gs_orig.variance);

    // residual state
    CHECK(ckpt.state.residual().variance == state.residual().variance);
    CHECK(ckpt.state.residual().y_adj.isApprox(state.residual().y_adj));

    // rng state: both rngs should produce the same next value
    const auto val_orig = rng();
    const auto val_read = ckpt.rng();
    CHECK(val_orig == val_read);

    std::filesystem::remove(ckpt_path);
}

TEST_CASE("checkpoint atomic write leaves no tmp file", "[checkpoint]")
{
    const std::string prefix = "/tmp/gelex_test_ckpt_atomic";
    const std::string ckpt_path = prefix + ".ckpt";
    const std::string tmp_path = prefix + ".ckpt.tmp";

    auto [model, priors] = make_bayes_a_model(kNSamples, kNSnps);
    BayesState state(model, priors);
    fill_state(state);
    std::mt19937_64 rng(99);

    write_checkpoint(state, rng, priors, prefix);

    REQUIRE(std::filesystem::exists(ckpt_path));
    REQUIRE_FALSE(std::filesystem::exists(tmp_path));

    std::filesystem::remove(ckpt_path);
}

TEST_CASE("checkpoint prior round-trip preserves all fields", "[checkpoint]")
{
    const std::string prefix = "/tmp/gelex_test_ckpt_prior";
    const std::string ckpt_path = prefix + ".ckpt";

    auto [model, model_priors] = make_bayes_a_model(kNSamples, kNSnps);
    std::mt19937_64 rng(77);

    Eigen::VectorXd pi_init(3);
    pi_init << 0.9, 0.05, 0.05;
    Eigen::VectorXd multiplier(3);
    multiplier << 0.0, 0.01, 0.1;

    std::vector<bayes::GeneticPrior> genetics;
    genetics.push_back(
        {.type = GeneticMode::A,
         .marker = bayes::
             MixturePrior{.variance = {.param = {.nu = 4.0, .s2 = 0.3}, .init = 0, .size = kNSnps}, .proportion = {.init = pi_init, .estimate = true}, .multiplier = multiplier},
         .sign = bayes::SignPrior{.init_value = 0.6}});

    std::vector<bayes::RandomPrior> random_priors;
    random_priors.push_back(
        {.param = {.nu = 3.0, .s2 = 0.1}, .init = 0, .size = 0});
    random_priors.push_back(
        {.param = {.nu = 5.0, .s2 = 0.2}, .init = 0, .size = 0});

    bayes::Priors priors(
        std::move(genetics),
        std::move(random_priors),
        {.param = {.nu = 4.0, .s2 = 0.5}, .init = 0, .size = 0});

    BayesState state(model, priors);
    fill_state(state);

    write_checkpoint(state, rng, priors, prefix);

    auto ckpt = read_checkpoint(ckpt_path);
    const auto& rp = ckpt.priors;

    // Residual
    CHECK(rp.residual().param.nu == 4.0);
    CHECK(rp.residual().param.s2 == 0.5);

    // Random
    REQUIRE(rp.random().size() == 2);
    CHECK(rp.random()[0].param.nu == 3.0);
    CHECK(rp.random()[0].param.s2 == 0.1);
    CHECK(rp.random()[1].param.nu == 5.0);
    CHECK(rp.random()[1].param.s2 == 0.2);

    // Genetic: MixturePrior
    REQUIRE(rp.genetics().size() == 1);
    CHECK(rp.genetics()[0].type == GeneticMode::A);
    const auto* mp = std::get_if<bayes::MixturePrior>(&rp.genetics()[0].marker);
    REQUIRE(mp != nullptr);
    CHECK(mp->variance.param.nu == 4.0);
    CHECK(mp->variance.param.s2 == 0.3);
    CHECK(mp->proportion.init.isApprox(pi_init));
    CHECK(mp->proportion.estimate == true);
    CHECK(mp->multiplier.isApprox(multiplier));

    // Sign
    REQUIRE(rp.genetics()[0].sign.has_value());
    CHECK(rp.genetics()[0].sign->init_value == 0.6);

    std::filesystem::remove(ckpt_path);
}

TEST_CASE("MCMC resume throws on checkpoint dimension mismatch", "[checkpoint]")
{
    const std::string prefix = "/tmp/gelex_test_ckpt_dimcheck";
    const std::string ckpt_path = prefix + ".ckpt";

    auto [model_write, priors_write] = make_bayes_a_model(kNSamples, kNSnps);
    BayesState state(model_write, priors_write);
    fill_state(state);
    std::mt19937_64 rng(42);

    write_checkpoint(state, rng, priors_write, prefix);

    auto ckpt = read_checkpoint(ckpt_path);

    auto [model_mismatch, priors_mismatch]
        = make_bayes_a_model(kNSamples, kNSnps + 3);
    (void)priors_mismatch;

    MCMCParams params(1, 0, 1, 1);
    MCMC mcmc(params, BayesA{});
    REQUIRE_THROWS(mcmc.resume(model_mismatch, std::move(ckpt)));

    std::filesystem::remove(ckpt_path);
}
