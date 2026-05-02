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
#include <cstdint>
#include <optional>
#include <random>
#include <span>
#include <utility>
#include <vector>

#include <Eigen/Core>
#include <catch2/catch_template_test_macros.hpp>
#include <catch2/catch_test_macros.hpp>

#include "gelex/algo/infer/mcmc/context.h"
#include "gelex/algo/infer/mcmc/kernels/bayes_a.h"
#include "gelex/algo/infer/mcmc/kernels/bayes_b.h"
#include "gelex/algo/infer/mcmc/kernels/bayes_c.h"
#include "gelex/algo/infer/mcmc/kernels/bayes_r.h"
#include "gelex/algo/infer/mcmc/kernels/bayes_rr.h"
#include "gelex/algo/infer/mcmc/recipes.h"
#include "gelex/algo/infer/mcmc/steps/fixed.h"
#include "gelex/algo/infer/mcmc/steps/genetic.h"
#include "gelex/algo/infer/mcmc/steps/pi.h"
#include "gelex/algo/infer/mcmc/steps/random.h"
#include "gelex/algo/infer/mcmc/sweep.h"
#include "gelex/data/genotype/genotype_matrix.h"
#include "gelex/exception.h"
#include "gelex/model/bayes/effects.h"
#include "gelex/model/bayes/genotype_storage.h"
#include "gelex/model/bayes/model.h"
#include "gelex/model/bayes/prior.h"
#include "gelex/model/bayes/states.h"
#include "gelex/types/fixed_effects.h"
#include "gelex/types/genetic_effect_type.h"

namespace gelex::mcmc
{
namespace
{

constexpr std::uint64_t kSeed = 0xC0FFEE5678ULL;

TEST_CASE(
    "FixedStep keeps residual identity and recovers OLS",
    "[mcmc][fixed]")
{
    const Eigen::MatrixXd X{
        {1, 0.5},
        {1, -0.3},
        {1, 0.1},
        {1, 0.7},
        {1, -0.4},
        {1, 0.2},
        {1, 0.9},
        {1, -0.6},
    };
    const Eigen::Vector2d beta{1.5, -2.0};
    const Eigen::VectorXd y = X * beta;

    FixedEffect effect;
    effect.X = X;
    effect.XtX_diag = X.colwise().squaredNorm();

    bayes::FixedState state{Eigen::VectorXd::Zero(X.cols())};
    bayes::ResidualState residual{.y_adj = y, .variance = 1e-6};
    std::mt19937_64 rng{kSeed};
    FixedStep sampler{FixedStep::Deps{
        .effect = effect,
        .state = state,
        .residual = residual,
        .rng = rng,
    }};

    for (int iter = 0; iter < 200; ++iter)
    {
        sampler.step();
    }

    REQUIRE((residual.y_adj + (X * state.coeffs)).isApprox(y));
    REQUIRE(state.coeffs.isApprox(beta, 1e-2));
}

TEST_CASE(
    "RandomStep keeps residual identity and updates variance",
    "[mcmc][random]")
{
    const Eigen::MatrixXd X{
        {1, 0, 0},
        {1, 0, 0},
        {0, 1, 0},
        {0, 1, 0},
        {0, 0, 1},
        {0, 0, 1},
    };
    const Eigen::Index n = X.rows();
    const Eigen::Index p = X.cols();

    const Eigen::VectorXd y = Eigen::VectorXd::LinSpaced(n, -1.0, 1.0);
    const bayes::RandomPrior prior{
        .param = {.nu = 4.0, .s2 = 0.5}, .init = 0.5, .size = 1};

    std::array<bayes::RandomEffect, 1> effects_arr{
        bayes::RandomEffect{"b", {"a", "b", "c"}, Eigen::MatrixXd{X}}};
    std::array<bayes::RandomState, 1> states_arr{
        bayes::RandomState{Eigen::VectorXd::Zero(p), prior.init}};
    bayes::ResidualState residual{.y_adj = y, .variance = 0.25};
    std::mt19937_64 rng{kSeed};

    RandomStep sampler{RandomStep::Deps{
        .effects = std::span<const bayes::RandomEffect>{effects_arr},
        .prior = prior,
        .states = std::span<bayes::RandomState>{states_arr},
        .residual = residual,
        .rng = rng,
    }};
    sampler.step();

    const auto& state = states_arr[0];
    REQUIRE(state.variance > 0.0);
    REQUIRE(std::isfinite(state.variance));
    REQUIRE((residual.y_adj + (X * state.coeffs)).isApprox(y));
}

constexpr Eigen::Index kTestP = 6;

template <typename Kernel>
struct KernelTraits;

template <>
struct KernelTraits<BayesAKernel>
{
    static auto make_prior() -> bayes::GeneticPrior
    {
        return bayes::GeneticPrior{
            .type = GeneticMode::A,
            .marker = bayes::
                ContinuousPrior{.variance = {.param = {.nu = 4.0, .s2 = 0.5}, .init = 0.1, .size = kTestP}},
            .sign = std::nullopt};
    }
};

template <>
struct KernelTraits<BayesRRKernel>
{
    static auto make_prior() -> bayes::GeneticPrior
    {
        return bayes::GeneticPrior{
            .type = GeneticMode::A,
            .marker = bayes::
                ContinuousPrior{.variance = {.param = {.nu = 4.0, .s2 = 0.5}, .init = 0.1, .size = 1}},
            .sign = std::nullopt};
    }
};

template <>
struct KernelTraits<BayesBKernel>
{
    static auto make_prior() -> bayes::GeneticPrior
    {
        return bayes::GeneticPrior{
            .type = GeneticMode::A,
            .marker = bayes::
                SpikePrior{.variance = {.param = {.nu = 4.0, .s2 = 0.5}, .init = 0.1, .size = kTestP}, .proportion = {.init = Eigen::VectorXd{{0.9, 0.1}}, .estimate = false}},
            .sign = std::nullopt};
    }
};

template <>
struct KernelTraits<BayesCKernel>
{
    static auto make_prior() -> bayes::GeneticPrior
    {
        return bayes::GeneticPrior{
            .type = GeneticMode::A,
            .marker = bayes::
                SpikePrior{.variance = {.param = {.nu = 4.0, .s2 = 0.5}, .init = 0.1, .size = 1}, .proportion = {.init = Eigen::VectorXd{{0.9, 0.1}}, .estimate = false}},
            .sign = std::nullopt};
    }
};

namespace
{

auto make_genetic_effect(Eigen::MatrixXd&& X) -> bayes::GeneticEffect
{
    const auto p = X.cols();
    Eigen::VectorXd mean = X.colwise().mean();
    Eigen::VectorXd stddev(p);
    for (Eigen::Index j = 0; j < p; ++j)
    {
        stddev(j) = std::sqrt(
            (X.col(j).array() - mean(j)).square().sum()
            / static_cast<double>(X.rows()));
    }
    gelex::genotype::GenotypeMatrix gm{
        std::move(X),
        std::vector<int64_t>{},
        std::move(mean),
        std::move(stddev)};
    return bayes::GeneticEffect{
        GeneticMode::A, bayes::GenotypeStorage{std::move(gm)}};
}

}  // namespace

TEMPLATE_TEST_CASE(
    "Genetic kernel preserves residual/GEBV identity and finite variance",
    "[mcmc][genetic]",
    BayesAKernel,
    BayesRRKernel,
    BayesBKernel,
    BayesCKernel)
{
    Eigen::MatrixXd X{
        {1, 0, 1, 0, 1, 0},
        {0, 1, 1, 0, 0, 1},
        {1, 1, 0, 1, 0, 0},
        {0, 0, 1, 1, 1, 1},
        {1, 0, 0, 1, 1, 0},
        {0, 1, 0, 0, 1, 1},
        {1, 1, 1, 0, 0, 0},
        {0, 0, 0, 1, 1, 1},
    };
    static_assert(kTestP == 6);
    const Eigen::VectorXd beta_true{{0.8, 0.0, -0.5, 0.0, 0.6, 0.0}};
    const Eigen::VectorXd y = X * beta_true;

    auto effect = make_genetic_effect(Eigen::MatrixXd{X});
    const auto prior = KernelTraits<TestType>::make_prior();
    bayes::GeneticState state{effect, prior};
    bayes::ResidualState residual{.y_adj = y, .variance = 1e-3};
    std::mt19937_64 rng{kSeed};

    GeneticStep<TestType> sampler{typename GeneticStep<TestType>::Deps{
        .block = {.effect = effect, .prior = prior, .state = state},
        .residual = residual,
        .rng = rng,
    }};

    for (int iter = 0; iter < 200; ++iter)
    {
        sampler.step();
    }

    REQUIRE((residual.y_adj + (X * state.coeffs)).isApprox(y));
    REQUIRE(state.u.isApprox(X * state.coeffs));
    REQUIRE(std::isfinite(state.variance));
    REQUIRE(state.variance >= 0.0);
    for (Eigen::Index i = 0; i < state.marker_variance.size(); ++i)
    {
        REQUIRE(std::isfinite(state.marker_variance(i)));
        REQUIRE(state.marker_variance(i) > 0.0);
    }
}

TEST_CASE(
    "Genetic kernel rejects mismatched prior/state",
    "[mcmc][genetic][errors]")
{
    Eigen::MatrixXd X{
        {1, 0, 1},
        {0, 1, 1},
        {1, 1, 0},
        {0, 0, 1},
    };
    const Eigen::Index p = X.cols();
    auto effect = make_genetic_effect(Eigen::MatrixXd{X});

    SECTION("BayesB rejects ContinuousPrior")
    {
        const auto prior = KernelTraits<BayesAKernel>::make_prior();
        bayes::GeneticState state{effect, prior};
        REQUIRE_THROWS_AS(BayesBKernel(prior, state), GelexException);
    }

    SECTION("BayesC rejects empty state.group")
    {
        const auto prior = KernelTraits<BayesCKernel>::make_prior();
        bayes::GeneticState state{
            GeneticMode::A,
            Eigen::VectorXd::Zero(p),
            Eigen::VectorXd::Zero(X.rows()),
            0.0,
            Eigen::VectorXd::Constant(1, 0.1),
            std::nullopt,
            std::nullopt};
        REQUIRE_THROWS_AS(BayesCKernel(prior, state), GelexException);
    }
}

// ──────────────────────────────────────────────────────────────────────────
// BayesR tests
// ──────────────────────────────────────────────────────────────────────────

namespace
{

auto make_mixture_prior(Eigen::Index p) -> bayes::GeneticPrior
{
    return bayes::
        GeneticPrior{
            .type = GeneticMode::A,
            .marker = bayes::
                MixturePrior{.variance = {.param = {.nu = 4.0, .s2 = 0.1}, .init = 0.05, .size = 1}, .proportion = {.init = Eigen::VectorXd{{0.7, 0.2, 0.08, 0.02}}, .estimate = true}, .multiplier = Eigen::VectorXd{{0.0, 0.001, 0.01, 0.1}}},
            .sign = std::nullopt};
    (void)p;
}

}  // namespace

TEST_CASE(
    "BayesRKernel sweep updates coeffs and component variance",
    "[mcmc][bayes-r]")
{
    constexpr Eigen::Index kN = 10;
    constexpr Eigen::Index kP = 8;
    constexpr Eigen::Index kK = 4;

    Eigen::MatrixXd X(kN, kP);
    X << 1, 0, 1, 0, 1, 0, 1, 0, 0, 1, 1, 0, 0, 1, 0, 1, 1, 1, 0, 1, 0, 0, 1, 0,
        0, 0, 1, 1, 1, 1, 0, 0, 1, 0, 0, 1, 1, 0, 0, 1, 0, 1, 0, 0, 1, 1, 0, 1,
        1, 1, 1, 0, 0, 0, 0, 1, 0, 0, 0, 1, 1, 1, 1, 0, 1, 0, 0, 1, 0, 1, 0, 0,
        0, 1, 0, 1, 1, 1, 0, 0;

    const Eigen::VectorXd beta_true{{0.8, 0.0, -0.5, 0.0, 0.6, 0.0, 0.0, 0.0}};
    const Eigen::VectorXd y = X * beta_true;

    auto effect = make_genetic_effect(Eigen::MatrixXd{X});
    const auto prior = make_mixture_prior(kP);
    bayes::GeneticState state{effect, prior};
    bayes::ResidualState residual{.y_adj = y, .variance = 1e-3};
    std::mt19937_64 rng{kSeed};

    GeneticSweep sweep{effect, state, residual, rng};
    BayesRKernel kernel{prior, state};
    sweep.run(kernel);

    auto& alloc = std::get<bayes::ComponentAllocation>(*state.group);

    REQUIRE(state.coeffs.size() == kP);
    REQUIRE(alloc.assignment.count.sum() == kP);
    for (Eigen::Index k = 0; k < kK; ++k)
    {
        REQUIRE(alloc.assignment.count(k) >= 0);
    }
    REQUIRE(state.marker_variance(0) > 0.0);
    REQUIRE(static_cast<Eigen::Index>(alloc.component_u.size()) == kK - 1);
    for (const auto& cu : alloc.component_u)
    {
        REQUIRE(cu.size() == kN);
    }
    REQUIRE(alloc.component_variance.size() == kK - 1);
}

TEST_CASE("BayesRKernel rejects non-MixturePrior", "[mcmc][bayes-r]")
{
    Eigen::MatrixXd X{
        {1, 0, 1},
        {0, 1, 1},
        {1, 1, 0},
        {0, 0, 1},
    };
    auto effect = make_genetic_effect(Eigen::MatrixXd{X});
    const auto prior = KernelTraits<BayesAKernel>::make_prior();
    bayes::GeneticState state{effect, prior};
    REQUIRE_THROWS_AS(BayesRKernel(prior, state), GelexException);
}

// ──────────────────────────────────────────────────────────────────────────
// PiStep tests
// ──────────────────────────────────────────────────────────────────────────

namespace
{

auto make_context_for_pi(
    const Eigen::MatrixXd& X,
    const Eigen::VectorXd& y,
    const bayes::GeneticPrior& genetic_prior,
    const bayes::ResidualPrior& residual_prior,
    mcmc::State& inference_state,
    std::mt19937_64& rng) -> std::pair<BayesModel, bayes::Priors>
{
    std::vector<bayes::GeneticEffect> genetics;
    genetics.emplace_back(
        make_genetic_effect(Eigen::MatrixXd{X}).type,
        std::move(make_genetic_effect(Eigen::MatrixXd{X}).X));

    BayesModel model{y, FixedEffect::build(y.size()), std::move(genetics)};

    bayes::Priors priors{
        std::vector<bayes::GeneticPrior>{genetic_prior},
        std::vector<bayes::RandomPrior>{},
        residual_prior};

    return {std::move(model), std::move(priors)};
}

}  // namespace

TEST_CASE(
    "PiStep updates proportion via Dirichlet posterior",
    "[mcmc][pi-sampler]")
{
    constexpr Eigen::Index kN = 10;
    constexpr Eigen::Index kP = 8;
    constexpr Eigen::Index kK = 4;

    Eigen::MatrixXd X(kN, kP);
    X << 1, 0, 1, 0, 1, 0, 1, 0, 0, 1, 1, 0, 0, 1, 0, 1, 1, 1, 0, 1, 0, 0, 1, 0,
        0, 0, 1, 1, 1, 1, 0, 0, 1, 0, 0, 1, 1, 0, 0, 1, 0, 1, 0, 0, 1, 1, 0, 1,
        1, 1, 1, 0, 0, 0, 0, 1, 0, 0, 0, 1, 1, 1, 1, 0, 1, 0, 0, 1, 0, 1, 0, 0,
        0, 1, 0, 1, 1, 1, 0, 0;

    const Eigen::VectorXd y = Eigen::VectorXd::LinSpaced(kN, -1.0, 1.0);
    const auto genetic_prior = make_mixture_prior(kP);

    auto effect = make_genetic_effect(Eigen::MatrixXd{X});
    bayes::GeneticState gstate{effect, genetic_prior};

    const bayes::ResidualPrior res_prior{
        .param = {.nu = 4.0, .s2 = 0.5}, .init = 0.5, .size = 1};

    std::vector<bayes::GeneticEffect> genetics_vec;
    genetics_vec.emplace_back(GeneticMode::A, std::move(effect.X));

    BayesModel model{y, FixedEffect::build(kN), std::move(genetics_vec)};

    bayes::Priors priors{
        std::vector<bayes::GeneticPrior>{genetic_prior},
        std::vector<bayes::RandomPrior>{},
        res_prior};

    mcmc::State inference_state{
        bayes::FixedState{Eigen::VectorXd::Zero(1)},
        std::vector<bayes::RandomState>{},
        std::vector<bayes::GeneticState>{std::move(gstate)},
        bayes::ResidualState{.y_adj = y, .variance = 0.5}};

    std::mt19937_64 rng{kSeed};
    Context ctx{.model = model, .priors = priors, .state = inference_state, .rng = rng};

    auto sampler = PiStep::make(ctx, GeneticMode::A);

    auto* gstate_ptr = inference_state.genetic(GeneticMode::A);
    auto& alloc = std::get<bayes::ComponentAllocation>(*gstate_ptr->group);
    alloc.assignment.count = Eigen::VectorXi{{6, 1, 1, 0}};

    sampler.step();

    const auto& prop = alloc.assignment.proportion;
    REQUIRE(prop.size() == kK);
    REQUIRE(std::abs(prop.sum() - 1.0) < 1e-10);
    for (Eigen::Index k = 0; k < kK; ++k)
    {
        REQUIRE(prop(k) > 0.0);
    }
    // component 0 (count=6) should dominate
    REQUIRE(prop(0) > prop(1));
    REQUIRE(prop(0) > prop(2));
}

TEST_CASE("PiStep rejects estimate=false prior", "[mcmc][pi-sampler]")
{
    constexpr Eigen::Index kN = 4;
    constexpr Eigen::Index kP = 3;

    Eigen::MatrixXd X{
        {1, 0, 1},
        {0, 1, 1},
        {1, 1, 0},
        {0, 0, 1},
    };
    const Eigen::VectorXd y = Eigen::VectorXd::Zero(kN);

    const bayes::GeneticPrior spike_prior{
        .type = GeneticMode::A,
        .marker = bayes::
            SpikePrior{.variance = {.param = {.nu = 4.0, .s2 = 0.5}, .init = 0.1, .size = 1}, .proportion = {.init = Eigen::VectorXd{{0.9, 0.1}}, .estimate = false}},
        .sign = std::nullopt};

    auto effect = make_genetic_effect(Eigen::MatrixXd{X});
    bayes::GeneticState gstate{effect, spike_prior};

    const bayes::ResidualPrior res_prior{
        .param = {.nu = 4.0, .s2 = 0.5}, .init = 0.5, .size = 1};

    std::vector<bayes::GeneticEffect> genetics_vec;
    genetics_vec.emplace_back(GeneticMode::A, std::move(effect.X));

    BayesModel model{y, FixedEffect::build(kN), std::move(genetics_vec)};

    bayes::Priors priors{
        std::vector<bayes::GeneticPrior>{spike_prior},
        std::vector<bayes::RandomPrior>{},
        res_prior};

    mcmc::State inference_state{
        bayes::FixedState{Eigen::VectorXd::Zero(1)},
        std::vector<bayes::RandomState>{},
        std::vector<bayes::GeneticState>{std::move(gstate)},
        bayes::ResidualState{.y_adj = y, .variance = 0.5}};

    std::mt19937_64 rng{kSeed};
    Context ctx{.model = model, .priors = priors, .state = inference_state, .rng = rng};

    REQUIRE_THROWS_AS(PiStep::make(ctx, GeneticMode::A), GelexException);
    (void)kP;
}

// ──────────────────────────────────────────────────────────────────────────
// BayesCpi / BayesBpi integration smoke tests
// ──────────────────────────────────────────────────────────────────────────

namespace
{

auto make_spike_prior_estimate(Eigen::Index size) -> bayes::GeneticPrior
{
    return bayes::GeneticPrior{
        .type = GeneticMode::A,
        .marker = bayes::
            SpikePrior{.variance = {.param = {.nu = 4.0, .s2 = 0.5}, .init = 0.1, .size = size}, .proportion = {.init = Eigen::VectorXd{{0.9, 0.1}}, .estimate = true}},
        .sign = std::nullopt};
}

auto make_chain_context(
    const Eigen::MatrixXd& X,
    const Eigen::VectorXd& y,
    const bayes::GeneticPrior& genetic_prior,
    mcmc::State& inference_state,
    std::mt19937_64& rng) -> std::pair<BayesModel, bayes::Priors>
{
    std::vector<bayes::GeneticEffect> genetics_vec;
    genetics_vec.emplace_back(
        GeneticMode::A, make_genetic_effect(Eigen::MatrixXd{X}).X);

    BayesModel model{y, FixedEffect::build(y.size()), std::move(genetics_vec)};

    const bayes::ResidualPrior res_prior{
        .param = {.nu = 4.0, .s2 = 0.5}, .init = 0.3, .size = 1};
    bayes::Priors priors{
        std::vector<bayes::GeneticPrior>{genetic_prior},
        std::vector<bayes::RandomPrior>{},
        res_prior};

    return {std::move(model), std::move(priors)};
}

auto build_inference_state(
    const bayes::GeneticPrior& genetic_prior,
    const Eigen::MatrixXd& X,
    const Eigen::VectorXd& y) -> mcmc::State
{
    auto effect = make_genetic_effect(Eigen::MatrixXd{X});
    bayes::GeneticState gstate{effect, genetic_prior};
    return mcmc::State{
        bayes::FixedState{Eigen::VectorXd::Zero(1)},
        std::vector<bayes::RandomState>{},
        std::vector<bayes::GeneticState>{std::move(gstate)},
        bayes::ResidualState{.y_adj = y, .variance = 0.3}};
}

}  // namespace

TEST_CASE("make_bayes_cpi_chain runs one step", "[mcmc][bayes-cpi]")
{
    Eigen::MatrixXd X{
        {1, 0, 1, 0, 1, 0},
        {0, 1, 1, 0, 0, 1},
        {1, 1, 0, 1, 0, 0},
        {0, 0, 1, 1, 1, 1},
        {1, 0, 0, 1, 1, 0},
        {0, 1, 0, 0, 1, 1},
        {1, 1, 1, 0, 0, 0},
        {0, 0, 0, 1, 1, 1},
    };
    const Eigen::VectorXd y
        = X * Eigen::VectorXd{{0.8, 0.0, -0.5, 0.0, 0.6, 0.0}};
    const auto genetic_prior = make_spike_prior_estimate(1);

    mcmc::State inference_state = build_inference_state(genetic_prior, X, y);
    std::mt19937_64 rng{kSeed};
    auto [model, priors]
        = make_chain_context(X, y, genetic_prior, inference_state, rng);

    Context ctx{.model = model, .priors = priors, .state = inference_state, .rng = rng};
    auto chain = make_bayes_cpi_chain(ctx);
    chain.step();

    const auto* gstate = inference_state.genetic(GeneticMode::A);
    REQUIRE(gstate != nullptr);
    REQUIRE(std::isfinite(gstate->variance));

    const auto& alloc = std::get<bayes::Assignment>(*gstate->group);
    REQUIRE(alloc.proportion.sum() > 0.99);
    REQUIRE(alloc.proportion.sum() < 1.01);
    REQUIRE(alloc.count.sum() == X.cols());
}

TEST_CASE("make_bayes_bpi_chain runs one step", "[mcmc][bayes-bpi]")
{
    Eigen::MatrixXd X{
        {1, 0, 1, 0, 1, 0},
        {0, 1, 1, 0, 0, 1},
        {1, 1, 0, 1, 0, 0},
        {0, 0, 1, 1, 1, 1},
        {1, 0, 0, 1, 1, 0},
        {0, 1, 0, 0, 1, 1},
        {1, 1, 1, 0, 0, 0},
        {0, 0, 0, 1, 1, 1},
    };
    const Eigen::VectorXd y
        = X * Eigen::VectorXd{{0.8, 0.0, -0.5, 0.0, 0.6, 0.0}};
    const auto genetic_prior = make_spike_prior_estimate(kTestP);

    mcmc::State inference_state = build_inference_state(genetic_prior, X, y);
    std::mt19937_64 rng{kSeed};
    auto [model, priors]
        = make_chain_context(X, y, genetic_prior, inference_state, rng);

    Context ctx{.model = model, .priors = priors, .state = inference_state, .rng = rng};
    auto chain = make_bayes_bpi_chain(ctx);
    chain.step();

    const auto* gstate = inference_state.genetic(GeneticMode::A);
    REQUIRE(gstate != nullptr);
    REQUIRE(std::isfinite(gstate->variance));

    const auto& alloc = std::get<bayes::Assignment>(*gstate->group);
    REQUIRE(alloc.proportion.sum() > 0.99);
    REQUIRE(alloc.proportion.sum() < 1.01);
    REQUIRE(alloc.count.sum() == X.cols());
}

}  // namespace
}  // namespace gelex::mcmc
