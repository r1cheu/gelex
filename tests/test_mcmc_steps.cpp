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
#include <memory>
#include <random>
#include <span>
#include <utility>
#include <vector>

#include <Eigen/Core>
#include <catch2/catch_template_test_macros.hpp>
#include <catch2/catch_test_macros.hpp>

#include "gelex/algo/infer/mcmc/context.h"
#include "gelex/algo/infer/mcmc/chain.h"
#include "gelex/algo/infer/mcmc/invariant.h"
#include "gelex/algo/infer/mcmc/kernels/bayes_a.h"
#include "gelex/algo/infer/mcmc/kernels/bayes_b.h"
#include "gelex/algo/infer/mcmc/kernels/bayes_c.h"
#include "gelex/algo/infer/mcmc/kernels/bayes_r.h"
#include "gelex/algo/infer/mcmc/kernels/bayes_rr.h"
#include "gelex/algo/infer/mcmc/recipes.h"
#include "gelex/algo/infer/mcmc/step.h"
#include "gelex/algo/infer/mcmc/steps/fixed.h"
#include "gelex/algo/infer/mcmc/steps/fixed_coefficient.h"
#include "gelex/algo/infer/mcmc/steps/genetic.h"
#include "gelex/algo/infer/mcmc/steps/pi.h"
#include "gelex/algo/infer/mcmc/steps/random.h"
#include "gelex/algo/infer/mcmc/steps/random_coefficient.h"
#include "gelex/algo/infer/mcmc/steps/random_variance.h"
#include "gelex/algo/infer/mcmc/steps/residual_variance.h"
#include "gelex/data/genotype/genotype.h"
#include "gelex/model/bayes/gaussian_prior.h"
#include "gelex/model/bayes/model.h"
#include "gelex/model/bayes/prior.h"
#include "gelex/model/bayes/prior_parameters.h"
#include "gelex/model/bayes/state.h"
#include "gelex/model/bayes/state_capabilities.h"
#include "gelex/types/fixed_designs.h"
#include "gelex/types/genetic_effect_type.h"
#include "genotype_fixture.h"

namespace gelex::mcmc
{
namespace
{

constexpr std::uint64_t kSeed = 0xC0FFEE5678ULL;

auto make_variance(double initial_value, double scale = 0.5)
    -> bayes::VarianceParameter
{
    return bayes::VarianceParameter(
        initial_value, bayes::ScaledInvChiSqPrior{4.0, scale});
}

auto make_random_prior(double initial_value) -> bayes::RandomPrior
{
    return bayes::RandomPrior{make_variance(initial_value)};
}

auto make_residual_prior(double initial_value) -> bayes::ResidualPrior
{
    return bayes::ResidualPrior{make_variance(initial_value)};
}

auto make_marker_variance(
    bayes::MarkerVarianceLayout scope,
    double initial_value = 0.1) -> bayes::MarkerVariance
{
    return bayes::MarkerVariance{scope, make_variance(initial_value)};
}

auto make_proportion(Eigen::VectorXd value, bayes::UpdatePolicy update)
    -> bayes::MixtureProportion
{
    const auto size = value.size();
    return bayes::MixtureProportion{
        bayes::SimplexParameter{
            std::move(value),
            bayes::DirichletPrior{Eigen::VectorXd::Ones(size)}},
        update};
}

auto make_genotype(Eigen::MatrixXd data) -> genotype::Genotype
{
    auto mean = data.colwise().mean().transpose().eval();
    auto stddev = Eigen::VectorXd::Ones(data.cols());
    return test::GenotypeBuilder::build(
        std::move(data), std::move(mean), std::move(stddev));
}

auto make_model(const Eigen::MatrixXd& X, const Eigen::VectorXd& y)
    -> BayesModel
{
    std::vector<bayes::GeneticDesign> genetics;
    genetics.emplace_back(GeneticMode::A, make_genotype(Eigen::MatrixXd{X}));
    return BayesModel{y, FixedDesign::build(y.size()), {}, std::move(genetics)};
}

template <typename Kernel>
auto make_prior() -> bayes::BayesPrior;

template <>
auto make_prior<BayesAKernel>() -> bayes::BayesPrior
{
    std::vector<std::unique_ptr<bayes::GeneticPrior>> genetics;
    genetics.push_back(
        std::make_unique<bayes::GaussianPrior>(
            GeneticMode::A,
            make_marker_variance(bayes::MarkerVarianceLayout::per_marker)));
    return bayes::BayesPrior(
        make_random_prior(0.5), std::move(genetics), make_residual_prior(0.25));
}

template <>
auto make_prior<BayesRRKernel>() -> bayes::BayesPrior
{
    std::vector<std::unique_ptr<bayes::GeneticPrior>> genetics;
    genetics.push_back(
        std::make_unique<bayes::GaussianPrior>(
            GeneticMode::A,
            make_marker_variance(bayes::MarkerVarianceLayout::shared)));
    return bayes::BayesPrior(
        make_random_prior(0.5), std::move(genetics), make_residual_prior(0.25));
}

template <>
auto make_prior<BayesBKernel>() -> bayes::BayesPrior
{
    std::vector<std::unique_ptr<bayes::GeneticPrior>> genetics;
    genetics.push_back(
        std::make_unique<bayes::SpikeSlabGaussianPrior>(
            GeneticMode::A,
            make_marker_variance(bayes::MarkerVarianceLayout::per_marker),
            make_proportion(
                Eigen::VectorXd{{0.9, 0.1}}, bayes::UpdatePolicy::fixed)));
    return bayes::BayesPrior(
        make_random_prior(0.5), std::move(genetics), make_residual_prior(0.25));
}

template <>
auto make_prior<BayesCKernel>() -> bayes::BayesPrior
{
    std::vector<std::unique_ptr<bayes::GeneticPrior>> genetics;
    genetics.push_back(
        std::make_unique<bayes::SpikeSlabGaussianPrior>(
            GeneticMode::A,
            make_marker_variance(bayes::MarkerVarianceLayout::shared),
            make_proportion(
                Eigen::VectorXd{{0.9, 0.1}}, bayes::UpdatePolicy::fixed)));
    return bayes::BayesPrior(
        make_random_prior(0.5), std::move(genetics), make_residual_prior(0.25));
}

auto make_bayes_r_prior() -> bayes::BayesPrior
{
    std::vector<std::unique_ptr<bayes::GeneticPrior>> genetics;
    genetics.push_back(
        std::make_unique<bayes::ScaledMixtureGaussianPrior>(
            GeneticMode::A,
            make_marker_variance(bayes::MarkerVarianceLayout::shared, 0.05),
            Eigen::VectorXd{{0.0, 0.001, 0.01, 0.1}},
            make_proportion(
                Eigen::VectorXd{{0.7, 0.2, 0.08, 0.02}},
                bayes::UpdatePolicy::sampled)));
    return bayes::BayesPrior(
        make_random_prior(0.5), std::move(genetics), make_residual_prior(0.25));
}

auto make_context(
    const BayesModel& model,
    const bayes::BayesPrior& prior,
    BayesState& state,
    std::mt19937_64& rng) -> Context
{
    return Context{.model = model, .prior = prior, .state = state, .rng = rng};
}

class RecordingStep final : public Step
{
   public:
    RecordingStep(std::vector<int>& calls, int value)
        : calls_(calls), value_(value)
    {
    }

    auto step() -> void override { calls_.push_back(value_); }

   private:
    std::vector<int>& calls_;
    int value_{};
};

}  // namespace

TEST_CASE("Runtime Chain runs heterogeneous steps in order", "[mcmc][chain]")
{
    std::vector<int> calls;
    std::vector<std::unique_ptr<Step>> steps;
    steps.push_back(std::make_unique<RecordingStep>(calls, 1));
    steps.push_back(std::make_unique<RecordingStep>(calls, 2));
    steps.push_back(std::make_unique<RecordingStep>(calls, 3));

    Chain chain{std::move(steps)};
    chain.step();

    const std::vector<int> expected{1, 2, 3};
    REQUIRE(calls == expected);
}

TEST_CASE(
    "ResidualAdjustmentGuard is no-op when coefficient is unchanged",
    "[mcmc][invariant]")
{
    const Eigen::VectorXd column{{1.0, 2.0, 3.0}};
    double coeff = 0.0;
    bayes::ResidualState residual{
        .y_adj = Eigen::VectorXd{{0.5, -0.25, 1.0}},
        .variance = 1.0};
    const auto before = residual.y_adj;

    {
        ResidualAdjustmentGuard guard{column, coeff, residual};
    }

    REQUIRE(residual.y_adj.isApprox(before));
}

TEST_CASE(
    "GeneticAdjustmentGuard updates residual and gebv in one transition",
    "[mcmc][invariant]")
{
    const Eigen::VectorXd column{{1.0, 2.0, -1.0}};
    double coeff = 0.5;
    bayes::ResidualState residual{
        .y_adj = Eigen::VectorXd{{1.0, 1.0, 1.0}},
        .variance = 1.0};
    bayes::GeneticState state{GeneticMode::A, 1, column.size()};

    {
        GeneticAdjustmentGuard guard{column, coeff, residual, state};
        coeff = -0.25;
    }

    const Eigen::VectorXd delta = (0.5 - coeff) * column;
    REQUIRE(residual.y_adj.isApprox(Eigen::VectorXd{{1.0, 1.0, 1.0}} + delta));
    REQUIRE(state.u.isApprox(-delta));
}

TEST_CASE(
    "GeneticMixtureAdjustmentGuard updates residual, gebv, and component gebv",
    "[mcmc][invariant]")
{
    const Eigen::VectorXd column{{1.0, 2.0, -1.0}};
    double coeff = 0.5;
    bayes::ResidualState residual{
        .y_adj = Eigen::VectorXd{{1.0, 1.0, 1.0}},
        .variance = 1.0};
    bayes::GeneticState state{GeneticMode::A, 1, column.size()};
    bayes::ComponentState component{2, column.size()};
    Eigen::VectorXi assignment{{1}};

    {
        GeneticMixtureAdjustmentGuard guard{
            column, coeff, residual, state, component, assignment, 0};
        coeff = -0.25;
        assignment(0) = 2;
    }

    const Eigen::VectorXd delta = (0.5 - coeff) * column;
    REQUIRE(residual.y_adj.isApprox(Eigen::VectorXd{{1.0, 1.0, 1.0}} + delta));
    REQUIRE(state.u.isApprox(-delta));
    REQUIRE(component.gebv[0].isApprox(-0.5 * column));
    REQUIRE(component.gebv[1].isApprox(coeff * column));
}

TEST_CASE(
    "FixedCoefficientStep keeps residual identity",
    "[mcmc][fixed][v2]")
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

    FixedDesign design;
    design.X = X;
    design.XtX_diag = X.colwise().squaredNorm();

    bayes::FixedState state{Eigen::VectorXd::Zero(X.cols())};
    bayes::ResidualState residual{.y_adj = y, .variance = 1e-6};
    std::mt19937_64 rng{kSeed};
    FixedCoefficientStep sampler{design, state, residual, rng};

    for (int iter = 0; iter < 200; ++iter)
    {
        sampler.step();
    }

    REQUIRE((residual.y_adj + (X * state.coeffs)).isApprox(y));
    REQUIRE(state.coeffs.isApprox(beta, 1e-2));
}

TEST_CASE(
    "Random coefficient and variance steps are independent",
    "[mcmc][random][v2]")
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

    std::array<bayes::RandomDesign, 1> designs{
        bayes::RandomDesign{"b", {"a", "b", "c"}, Eigen::MatrixXd{X}}};
    std::array<bayes::RandomState, 1> states{
        bayes::RandomState{Eigen::VectorXd::Zero(p), 0.5}};
    bayes::ResidualState residual{.y_adj = y, .variance = 0.25};
    auto prior = make_random_prior(0.5);
    std::mt19937_64 rng{kSeed};

    RandomCoefficientStep coeff_step{
        std::span<const bayes::RandomDesign>{designs},
        std::span<bayes::RandomState>{states},
        residual,
        rng};
    coeff_step.step();

    const double variance_before = states[0].variance;
    REQUIRE(states[0].coeffs.squaredNorm() > 0.0);
    REQUIRE(variance_before == 0.5);
    REQUIRE((residual.y_adj + (X * states[0].coeffs)).isApprox(y));

    const auto coeffs_before = states[0].coeffs;
    RandomVarianceStep variance_step{
        prior, std::span<bayes::RandomState>{states}, rng};
    variance_step.step();

    REQUIRE(states[0].coeffs.isApprox(coeffs_before));
    REQUIRE(states[0].variance > 0.0);
    REQUIRE(std::isfinite(states[0].variance));
}

TEST_CASE(
    "ResidualVarianceStep only updates residual variance",
    "[mcmc][residual][v2]")
{
    bayes::ResidualState residual{
        .y_adj = Eigen::VectorXd{{1.0, -0.5, 0.25, -0.75}},
        .variance = 0.25};
    const auto y_adj_before = residual.y_adj;
    auto prior = make_residual_prior(0.25);
    std::mt19937_64 rng{kSeed};

    ResidualVarianceStep sampler{residual.y_adj.size(), prior, residual, rng};
    sampler.step();

    REQUIRE(residual.y_adj.isApprox(y_adj_before));
    REQUIRE(residual.variance > 0.0);
    REQUIRE(std::isfinite(residual.variance));
}

TEST_CASE("FixedStep keeps residual identity and recovers OLS", "[mcmc][fixed]")
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

    FixedDesign design;
    design.X = X;
    design.XtX_diag = X.colwise().squaredNorm();

    bayes::FixedState state{Eigen::VectorXd::Zero(X.cols())};
    bayes::ResidualState residual{.y_adj = y, .variance = 1e-6};
    std::mt19937_64 rng{kSeed};
    FixedStep sampler{FixedStep::Deps{
        .design = design,
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

    std::array<bayes::RandomDesign, 1> designs{
        bayes::RandomDesign{"b", {"a", "b", "c"}, Eigen::MatrixXd{X}}};
    std::array<bayes::RandomState, 1> states{
        bayes::RandomState{Eigen::VectorXd::Zero(p), 0.5}};
    bayes::ResidualState residual{.y_adj = y, .variance = 0.25};
    std::mt19937_64 rng{kSeed};

    RandomStep sampler{RandomStep::Deps{
        .designs = std::span<const bayes::RandomDesign>{designs},
        .variance = make_random_prior(0.5),
        .states = std::span<bayes::RandomState>{states},
        .residual = residual,
        .rng = rng,
    }};
    sampler.step();

    REQUIRE(states[0].variance > 0.0);
    REQUIRE(std::isfinite(states[0].variance));
    REQUIRE((residual.y_adj + (X * states[0].coeffs)).isApprox(y));
}

TEMPLATE_TEST_CASE(
    "Independent genetic kernel preserves residual/GEBV identity",
    "[mcmc][genetic]",
    BayesAKernel,
    BayesRRKernel,
    BayesBKernel,
    BayesCKernel)
{
    const Eigen::MatrixXd X{
        {1, 0, 1, 0, 1, 0},
        {0, 1, 1, 0, 0, 1},
        {1, 1, 0, 1, 0, 0},
        {0, 0, 1, 1, 1, 1},
        {1, 0, 0, 1, 1, 0},
        {0, 1, 0, 0, 1, 1},
        {1, 1, 1, 0, 0, 0},
        {0, 0, 0, 1, 1, 1},
    };
    const Eigen::VectorXd beta{{0.8, 0.0, -0.5, 0.0, 0.6, 0.0}};
    const Eigen::VectorXd y = X * beta;

    auto model = make_model(X, y);
    auto prior = make_prior<TestType>();
    BayesState state(model, prior);
    std::mt19937_64 rng{kSeed};
    auto ctx = make_context(model, prior, state, rng);
    auto sampler = GeneticStep<TestType>::make(ctx, GeneticMode::A);

    for (int iter = 0; iter < 200; ++iter)
    {
        sampler.step();
    }

    const auto& genetic = *state.genetic(GeneticMode::A);
    REQUIRE((state.residual().y_adj + (X * genetic.coeffs)).isApprox(y));
    REQUIRE(genetic.u.isApprox(X * genetic.coeffs));
    REQUIRE(std::isfinite(genetic.variance));
    REQUIRE(genetic.variance >= 0.0);

    const auto& variance = state.genetic_block_for(GeneticMode::A)
                               ->prior_state()
                               .require<bayes::VarianceStateCap>()
                               .variance()[0];
    for (Eigen::Index i = 0; i < variance.size(); ++i)
    {
        REQUIRE(std::isfinite(variance(i)));
        REQUIRE(variance(i) > 0.0);
    }
}

TEST_CASE("BayesRKernel updates component state", "[mcmc][bayes-r]")
{
    const Eigen::MatrixXd X{
        {1, 0, 1, 0, 1, 0, 1, 0},
        {0, 1, 1, 0, 0, 1, 0, 1},
        {1, 1, 0, 1, 0, 0, 1, 0},
        {0, 0, 1, 1, 1, 1, 0, 0},
        {1, 0, 0, 1, 1, 0, 0, 1},
        {0, 1, 0, 0, 1, 1, 0, 1},
        {1, 1, 1, 0, 0, 0, 0, 1},
        {0, 0, 0, 1, 1, 1, 1, 0},
    };
    const Eigen::VectorXd beta{{0.8, 0.0, -0.5, 0.0, 0.6, 0.0, 0.0, 0.0}};
    const Eigen::VectorXd y = X * beta;

    auto model = make_model(X, y);
    auto prior = make_bayes_r_prior();
    BayesState state(model, prior);
    std::mt19937_64 rng{kSeed};
    auto ctx = make_context(model, prior, state, rng);
    auto sampler = GeneticStep<BayesRKernel>::make(ctx, GeneticMode::A);
    sampler.step();

    const auto& prior_state
        = state.genetic_block_for(GeneticMode::A)->prior_state();
    const auto& component
        = prior_state.require<bayes::ComponentStateCap>().component()[0];
    const auto& proportion
        = prior_state.require<bayes::ProportionStateCap>().proportion()[0];

    REQUIRE(proportion.count.sum() == X.cols());
    REQUIRE(
        static_cast<Eigen::Index>(component.gebv.size())
        == proportion.value.size() - 1);
    REQUIRE(component.gebv_var.size() == proportion.value.size() - 1);
}

TEST_CASE("PiStep samples or skips proportion by state update", "[mcmc][pi]")
{
    const Eigen::MatrixXd X{
        {1, 0, 1, 0},
        {0, 1, 1, 0},
        {1, 1, 0, 1},
        {0, 0, 1, 1},
    };
    const Eigen::VectorXd y = Eigen::VectorXd::LinSpaced(X.rows(), -1.0, 1.0);
    auto model = make_model(X, y);

    SECTION("sampled proportion is updated")
    {
        auto prior = make_bayes_r_prior();
        BayesState state(model, prior);
        auto& proportion = state.genetic_block_for(GeneticMode::A)
                               ->prior_state()
                               .require<bayes::ProportionStateCap>()
                               .proportion()[0];
        proportion.count = Eigen::VectorXi{{6, 1, 1, 0}};
        const auto before = proportion.value;

        std::mt19937_64 rng{kSeed};
        auto ctx = make_context(model, prior, state, rng);
        auto sampler = PiStep::make(ctx, GeneticMode::A);
        sampler.step();

        REQUIRE(proportion.value.size() == before.size());
        REQUIRE(std::abs(proportion.value.sum() - 1.0) < 1e-10);
        REQUIRE_FALSE(proportion.value.isApprox(before, 0.0));
    }

    SECTION("fixed proportion is no-op")
    {
        auto prior = make_prior<BayesCKernel>();
        BayesState state(model, prior);
        auto& proportion = state.genetic_block_for(GeneticMode::A)
                               ->prior_state()
                               .require<bayes::ProportionStateCap>()
                               .proportion()[0];
        proportion.count = Eigen::VectorXi{{2, 2}};
        const auto before = proportion.value;

        std::mt19937_64 rng{kSeed};
        auto ctx = make_context(model, prior, state, rng);
        auto sampler = PiStep::make(ctx, GeneticMode::A);
        sampler.step();

        REQUIRE(proportion.value.isApprox(before, 0.0));
    }
}

TEST_CASE("BayesCpi and BayesBpi chains run one step", "[mcmc][chain]")
{
    const Eigen::MatrixXd X{
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

    auto model = make_model(X, y);

    SECTION("BayesCpi")
    {
        auto prior = make_prior<BayesCKernel>();
        BayesState state(model, prior);
        std::mt19937_64 rng{kSeed};
        auto ctx = make_context(model, prior, state, rng);
        auto chain = make_bayes_cpi_chain<GeneticShape::a_only>(ctx);
        chain.step();
        REQUIRE(std::isfinite(state.genetic(GeneticMode::A)->variance));
    }

    SECTION("BayesBpi")
    {
        auto prior = make_prior<BayesBKernel>();
        BayesState state(model, prior);
        std::mt19937_64 rng{kSeed};
        auto ctx = make_context(model, prior, state, rng);
        auto chain = make_bayes_bpi_chain<GeneticShape::a_only>(ctx);
        chain.step();
        REQUIRE(std::isfinite(state.genetic(GeneticMode::A)->variance));
    }
}

}  // namespace gelex::mcmc
