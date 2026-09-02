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
#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>
#include <concepts>
#include <random>
#include <string>
#include <type_traits>
#include <utility>
#include <vector>

#include "gelex/algo/mcmc/chain.h"
#include "gelex/bayes/genetic/gaussian_prior.h"
#include "gelex/bayes/genetic/gaussian_prior_state.h"
#include "gelex/bayes/genetic/half_normal_prior.h"
#include "gelex/bayes/genetic/independent_topology.h"
#include "gelex/bayes/genetic/joint_topology.h"
#include "gelex/bayes/genetic/legacy_genetic_prior.h"
#include "gelex/bayes/genetic/parameters.h"
#include "gelex/bayes/kernel.h"
#include "gelex/bayes/legacy_prior.h"
#include "gelex/bayes/legacy_state.h"
#include "gelex/bayes/parameter/distributions.h"
#include "gelex/bayes/parameter/values.h"
#include "gelex/bayes/prior_compilation.h"
#include "gelex/bayes/recipe.h"
#include "gelex/bayes/semantic_method.h"
#include "gelex/bayes/spec.h"
#include "gelex/bayes/state.h"
#include "gelex/bayes/variance_budget.h"
#include "gelex/data/genotype_method.h"
#include "gelex/exception.h"
#include "gelex/types/fixed_designs.h"
#include "gelex/types/genetic_mode.h"

#include "compact_genotype_fixture.h"

using Catch::Approx;

namespace
{

constexpr auto mode_a = gelex::GeneticModeSet{gelex::GeneticMode::A};
constexpr auto mode_ad = gelex::GeneticMode::A | gelex::GeneticMode::D;
using Method = gelex::GaussianMethod<gelex::VarianceLayout::Pooled>;
using UnpooledMethod = gelex::GaussianMethod<gelex::VarianceLayout::Unpooled>;
using PooledSpikeSlabMethod
    = gelex::SpikeSlabMethod<gelex::VarianceLayout::Pooled>;
using UnpooledSpikeSlabMethod
    = gelex::SpikeSlabMethod<gelex::VarianceLayout::Unpooled>;
using FixedUnpooledSpikeSlabMethod = gelex::SpikeSlabMethod<
    gelex::VarianceLayout::Unpooled,
    gelex::UpdatePolicy::Fixed>;
using SampledScaledMixtureMethod = gelex::ScaledMixtureMethod<>;
using FixedScaledMixtureMethod
    = gelex::ScaledMixtureMethod<gelex::UpdatePolicy::Fixed>;
using SampledJointSpikeSlabMethod = gelex::JointSpikeSlabMethod<>;
using FixedJointSpikeSlabMethod = gelex::JointSpikeSlabMethod<
    gelex::UpdatePolicy::Fixed,
    gelex::UpdatePolicy::Fixed>;
using AdditiveGeneticPrior = gelex::IndependentTopology<
    mode_a,
    gelex::GaussianPrior<gelex::VarianceLayout::Pooled>>;
using AdditiveDominanceGeneticPrior = gelex::IndependentTopology<
    mode_ad,
    gelex::GaussianPrior<gelex::VarianceLayout::Pooled>>;
using UnpooledAdditiveGeneticPrior = gelex::IndependentTopology<
    mode_a,
    gelex::GaussianPrior<gelex::VarianceLayout::Unpooled>>;
using PooledSpikeSlabGeneticPrior = gelex::IndependentTopology<
    mode_a,
    gelex::SpikeSlabPrior<gelex::VarianceLayout::Pooled>>;
using PooledSpikeSlabADGeneticPrior = gelex::IndependentTopology<
    mode_ad,
    gelex::SpikeSlabPrior<gelex::VarianceLayout::Pooled>>;
using UnpooledSpikeSlabGeneticPrior = gelex::IndependentTopology<
    mode_a,
    gelex::SpikeSlabPrior<gelex::VarianceLayout::Unpooled>>;
using SampledScaledMixtureGeneticPrior
    = gelex::IndependentTopology<mode_a, gelex::ScaledMixturePrior<>>;
using FixedScaledMixtureGeneticPrior = gelex::IndependentTopology<
    mode_a,
    gelex::ScaledMixturePrior<gelex::UpdatePolicy::Fixed>>;
using SampledJointSpikeSlabGeneticPrior = gelex::JointTopology<
    gelex::GaussianPrior<gelex::VarianceLayout::Pooled>,
    gelex::JointSpikeSlabPrior<>>;
using FixedJointSpikeSlabGeneticPrior = gelex::JointTopology<
    gelex::GaussianPrior<gelex::VarianceLayout::Pooled>,
    gelex::JointSpikeSlabPrior<
        gelex::UpdatePolicy::Fixed,
        gelex::UpdatePolicy::Fixed>>;
using FixedAllocationJointSpikeSlabGeneticPrior = gelex::JointTopology<
    gelex::GaussianPrior<gelex::VarianceLayout::Pooled>,
    gelex::JointSpikeSlabPrior<
        gelex::UpdatePolicy::Fixed,
        gelex::UpdatePolicy::Sampled>>;
using FixedSignJointSpikeSlabGeneticPrior = gelex::JointTopology<
    gelex::GaussianPrior<gelex::VarianceLayout::Pooled>,
    gelex::JointSpikeSlabPrior<
        gelex::UpdatePolicy::Sampled,
        gelex::UpdatePolicy::Fixed>>;

static_assert(std::is_empty_v<
              gelex::detail::ProbabilityUpdater<gelex::UpdatePolicy::Fixed>>);
static_assert(!std::is_empty_v<
              gelex::detail::ProbabilityUpdater<gelex::UpdatePolicy::Sampled>>);
static_assert(std::is_empty_v<gelex::detail::SimplexUpdater<
                  gelex::UpdatePolicy::Fixed,
                  gelex::ScaledMixturePrior<>::class_count>>);
static_assert(!std::is_empty_v<gelex::detail::SimplexUpdater<
                  gelex::UpdatePolicy::Sampled,
                  gelex::ScaledMixturePrior<>::class_count>>);

static_assert(
    std::same_as<
        decltype(gelex::make_kernel(
            std::declval<const gelex::BayesPrior<AdditiveGeneticPrior>&>())),
        gelex::BayesKernel<AdditiveGeneticPrior>>);

static_assert(
    std::same_as<
        decltype(gelex::make_kernel(
            std::declval<
                const gelex::BayesPrior<AdditiveDominanceGeneticPrior>&>())),
        gelex::BayesKernel<AdditiveDominanceGeneticPrior>>);

static_assert(
    std::same_as<
        decltype(gelex::make_kernel(
            std::declval<
                const gelex::BayesPrior<UnpooledAdditiveGeneticPrior>&>())),
        gelex::BayesKernel<UnpooledAdditiveGeneticPrior>>);

static_assert(
    std::same_as<
        decltype(gelex::make_kernel(
            std::declval<
                const gelex::BayesPrior<PooledSpikeSlabGeneticPrior>&>())),
        gelex::BayesKernel<PooledSpikeSlabGeneticPrior>>);

static_assert(
    std::same_as<
        decltype(gelex::make_kernel(
            std::declval<
                const gelex::BayesPrior<PooledSpikeSlabADGeneticPrior>&>())),
        gelex::BayesKernel<PooledSpikeSlabADGeneticPrior>>);

static_assert(
    std::same_as<
        decltype(gelex::make_kernel(
            std::declval<
                const gelex::BayesPrior<UnpooledSpikeSlabGeneticPrior>&>())),
        gelex::BayesKernel<UnpooledSpikeSlabGeneticPrior>>);

static_assert(
    std::same_as<
        decltype(gelex::make_kernel(
            std::declval<
                const gelex::BayesPrior<SampledScaledMixtureGeneticPrior>&>())),
        gelex::BayesKernel<SampledScaledMixtureGeneticPrior>>);

static_assert(
    std::same_as<
        decltype(gelex::make_kernel(
            std::declval<
                const gelex::BayesPrior<FixedScaledMixtureGeneticPrior>&>())),
        gelex::BayesKernel<FixedScaledMixtureGeneticPrior>>);

static_assert(std::same_as<
              decltype(gelex::make_kernel(
                  std::declval<const gelex::BayesPrior<
                      SampledJointSpikeSlabGeneticPrior>&>())),
              gelex::BayesKernel<SampledJointSpikeSlabGeneticPrior>>);

static_assert(
    std::same_as<
        decltype(gelex::make_kernel(
            std::declval<
                const gelex::BayesPrior<FixedJointSpikeSlabGeneticPrior>&>())),
        gelex::BayesKernel<FixedJointSpikeSlabGeneticPrior>>);

static_assert(std::same_as<
              decltype(gelex::make_kernel(
                  std::declval<const gelex::BayesPrior<
                      FixedAllocationJointSpikeSlabGeneticPrior>&>())),
              gelex::BayesKernel<FixedAllocationJointSpikeSlabGeneticPrior>>);

static_assert(std::same_as<
              decltype(gelex::make_kernel(
                  std::declval<const gelex::BayesPrior<
                      FixedSignJointSpikeSlabGeneticPrior>&>())),
              gelex::BayesKernel<FixedSignJointSpikeSlabGeneticPrior>>);

auto make_model() -> gelex::BayesModel
{
    return gelex::test::make_compact_model(
        Eigen::MatrixXd{{0.0, 1.0}, {1.0, 1.0}, {2.0, 1.0}, {0.0, 1.0}},
        Eigen::VectorXd{{1.0, -0.5, 0.25, 2.0}},
        mode_a);
}

auto make_ad_model() -> gelex::BayesModel
{
    return gelex::test::make_compact_model(
        Eigen::MatrixXd{{0.0, 1.0}, {1.0, 1.0}, {2.0, 1.0}, {0.0, 1.0}},
        Eigen::VectorXd{{1.0, -0.5, 0.25, 2.0}},
        mode_ad);
}

auto make_joint_model() -> gelex::BayesModel
{
    return gelex::test::make_compact_model(
        Eigen::MatrixXd{{0.0, 1.0}, {1.0, 1.0}, {2.0, 1.0}, {0.0, 1.0}},
        Eigen::VectorXd{{1.0, -0.5, 0.25, 2.0}},
        mode_ad,
        gelex::GenotypeMethod::NOIACenter);
}

auto make_model_with_random() -> gelex::BayesModel
{
    auto genetic = gelex::test::make_genetic_design(
        Eigen::MatrixXd{{0.0, 1.0}, {1.0, 1.0}, {2.0, 1.0}, {0.0, 1.0}},
        mode_a);
    std::vector<gelex::bayes::RandomDesign> random;
    random.emplace_back(
        "batch",
        std::vector<std::string>{"batch"},
        Eigen::MatrixXd{{0.0}, {1.0}, {0.0}, {1.0}});
    return gelex::BayesModel{
        Eigen::VectorXd{{1.0, -0.5, 0.25, 2.0}},
        gelex::FixedDesign::make(4),
        std::move(random),
        std::move(genetic)};
}

auto reconstruct_genetic_fitted(
    const gelex::bayes::GeneticDesign& design,
    gelex::GeneticMode mode,
    const Eigen::VectorXd& coefficients) -> Eigen::VectorXd
{
    Eigen::VectorXd fitted = Eigen::VectorXd::Zero(design.rows());
    const auto& projection = design.projection(mode);
    for (Eigen::Index marker = 0; marker < coefficients.size(); ++marker)
    {
        projection.axpy(marker, coefficients(marker), fitted);
    }
    return fitted;
}

template <gelex::ComponentLayout Layout>
auto reconstruct_component_fitted(
    const gelex::bayes::GeneticDesign& design,
    gelex::GeneticMode mode,
    const Eigen::VectorXd& coefficients,
    const Eigen::VectorX<std::uint8_t>& assignment) -> Eigen::MatrixXd
{
    Eigen::MatrixXd fitted
        = Eigen::MatrixXd::Zero(design.rows(), Layout::component_count);
    const auto& projection = design.projection(mode);
    for (Eigen::Index marker = 0; marker < coefficients.size(); ++marker)
    {
        const int component = Layout::component_index(
            mode, static_cast<std::size_t>(assignment(marker)));
        if (component != Layout::no_component)
        {
            projection.axpy(
                marker, coefficients(marker), fitted.col(component));
        }
    }
    return fitted;
}

auto make_legacy_gaussian_prior(
    gelex::GeneticMode mode,
    const gelex::VarianceParameter& variance) -> gelex::bayes::GeneticPrior
{
    gelex::bayes::SingleGeneticPrior genetic{
        gelex::bayes::SingleSharedGaussianPrior{
            mode,
            gelex::bayes::SharedMarkerVariance{gelex::bayes::VarianceParameter{
                variance.initial_value(),
                gelex::bayes::ScaledInvChiSqPrior{
                    variance.prior().degrees_of_freedom(),
                    variance.prior().scale()}}}}};
    return gelex::bayes::GeneticPrior{std::move(genetic)};
}

auto make_legacy_unpooled_gaussian_prior(
    gelex::GeneticMode mode,
    const gelex::VarianceParameter& variance) -> gelex::bayes::GeneticPrior
{
    gelex::bayes::SingleGeneticPrior genetic{
        gelex::bayes::SinglePerMarkerGaussianPrior{
            mode,
            gelex::bayes::PerMarkerVariance{gelex::bayes::VarianceParameter{
                variance.initial_value(),
                gelex::bayes::ScaledInvChiSqPrior{
                    variance.prior().degrees_of_freedom(),
                    variance.prior().scale()}}}}};
    return gelex::bayes::GeneticPrior{std::move(genetic)};
}

auto make_legacy_joint_prior(const SampledJointSpikeSlabGeneticPrior& prior)
    -> gelex::bayes::GeneticPrior
{
    const auto make_variance = [](const gelex::VarianceParameter& parameter)
    {
        return gelex::bayes::SharedMarkerVariance{
            gelex::bayes::VarianceParameter{
                parameter.initial_value(),
                gelex::bayes::ScaledInvChiSqPrior{
                    parameter.prior().degrees_of_freedom(),
                    parameter.prior().scale()}}};
    };
    const auto& joint = prior.joint();
    const Eigen::VectorXd probabilities
        = Eigen::Map<const Eigen::Vector<double, 4>>{
            joint.probabilities.initial.data()};
    gelex::bayes::JointGeneticPrior genetic{
        gelex::bayes::JointHalfNormalMixturePrior{
            gelex::bayes::JointSharedMarkerVariance{std::array{
                make_variance(
                    prior.mode_values().get<gelex::GeneticMode::A>().variance),
                make_variance(prior.mode_values()
                                  .get<gelex::GeneticMode::D>()
                                  .variance)}},
            gelex::bayes::MixtureProportion{gelex::bayes::SimplexParameter{
                probabilities,
                gelex::bayes::DirichletPrior{Eigen::VectorXd::Ones(4)}}},
            gelex::bayes::ProbabilityParameter{
                joint.positive_probability.initial,
                gelex::bayes::BetaPrior{1.0, 1.0}}}};
    return gelex::bayes::GeneticPrior{std::move(genetic)};
}

template <typename GeneticPrior>
auto make_legacy_prior(
    const gelex::BayesPrior<GeneticPrior>& prior,
    std::vector<gelex::bayes::GeneticPrior> genetics)
    -> gelex::bayes::LegacyBayesPrior
{
    return gelex::bayes::LegacyBayesPrior{
        gelex::bayes::RandomPrior{
            1.0, gelex::bayes::ScaledInvChiSqPrior{4.0, 0.5}},
        std::move(genetics),
        gelex::bayes::ResidualPrior{
            prior.residual().initial_value(),
            gelex::bayes::ScaledInvChiSqPrior{
                prior.residual().prior().degrees_of_freedom(),
                prior.residual().prior().scale()}}};
}

}  // namespace

TEST_CASE(
    "pooled Gaussian kernel preserves adjusted-response and fitted-cache "
    "invariants",
    "[bayes][kernel]")
{
    const auto model = make_model();
    const auto prior
        = gelex::compile(gelex::BayesRecipe<Method, mode_a>::defaults(), model);
    auto state = gelex::make_state(prior, model);
    auto kernel = gelex::make_kernel(prior);
    std::mt19937_64 rng{123};

    kernel.step(model, state, rng);

    const auto& genetic = state.genetic().get<gelex::GeneticMode::A>();
    const auto reconstructed = reconstruct_genetic_fitted(
        model.genetic(), gelex::GeneticMode::A, genetic.coefficients);
    const Eigen::VectorXd fixed_fitted
        = model.fixed().X * state.fixed().coefficients;

    REQUIRE(genetic.component_fitted_values.col(0).isApprox(reconstructed));
    REQUIRE((state.residual().adjusted_response + fixed_fitted + reconstructed)
                .isApprox(model.phenotype()));
    REQUIRE(genetic.coefficients(1) == 0.0);
    REQUIRE(state.residual().variance > 0.0);
    REQUIRE(genetic.family_state.variance > 0.0);
}

TEST_CASE(
    "independent AD pooled Gaussian kernel preserves mode caches and the "
    "adjusted-response invariant",
    "[bayes][kernel]")
{
    const auto model = make_ad_model();
    const auto prior = gelex::compile(
        gelex::BayesRecipe<Method, mode_ad>::defaults(), model);
    auto state = gelex::make_state(prior, model);
    auto kernel = gelex::make_kernel(prior);
    std::mt19937_64 rng{123};

    kernel.step(model, state, rng);

    const auto& additive = state.genetic().get<gelex::GeneticMode::A>();
    const auto& dominance = state.genetic().get<gelex::GeneticMode::D>();
    const auto additive_fitted = reconstruct_genetic_fitted(
        model.genetic(), gelex::GeneticMode::A, additive.coefficients);
    const auto dominance_fitted = reconstruct_genetic_fitted(
        model.genetic(), gelex::GeneticMode::D, dominance.coefficients);
    const Eigen::VectorXd fixed_fitted
        = model.fixed().X * state.fixed().coefficients;

    REQUIRE(additive.component_fitted_values.col(0).isApprox(additive_fitted));
    REQUIRE(
        dominance.component_fitted_values.col(0).isApprox(dominance_fitted));
    REQUIRE((state.residual().adjusted_response + fixed_fitted + additive_fitted
             + dominance_fitted)
                .isApprox(model.phenotype()));
    REQUIRE(additive.coefficients(1) == 0.0);
    REQUIRE(dominance.coefficients(1) == 0.0);
    REQUIRE(additive.family_state.variance > 0.0);
    REQUIRE(dominance.family_state.variance > 0.0);
}

TEST_CASE(
    "unpooled Gaussian kernel preserves adjusted-response and fitted-cache "
    "invariants",
    "[bayes][kernel]")
{
    const auto model = make_model();
    const auto prior = gelex::compile(
        gelex::BayesRecipe<UnpooledMethod, mode_a>::defaults(), model);
    auto state = gelex::make_state(prior, model);
    auto kernel = gelex::make_kernel(prior);
    const double invalid_marker_variance
        = state.genetic().get<gelex::GeneticMode::A>().family_state.variance(1);
    std::mt19937_64 rng{123};

    kernel.step(model, state, rng);

    const auto& genetic = state.genetic().get<gelex::GeneticMode::A>();
    const auto reconstructed = reconstruct_genetic_fitted(
        model.genetic(), gelex::GeneticMode::A, genetic.coefficients);
    const Eigen::VectorXd fixed_fitted
        = model.fixed().X * state.fixed().coefficients;

    REQUIRE(genetic.component_fitted_values.col(0).isApprox(reconstructed));
    REQUIRE((state.residual().adjusted_response + fixed_fitted + reconstructed)
                .isApprox(model.phenotype()));
    REQUIRE(genetic.coefficients(1) == 0.0);
    REQUIRE(genetic.family_state.variance.size() == model.genetic().cols());
    REQUIRE((genetic.family_state.variance.array() > 0.0).all());
    REQUIRE(genetic.family_state.variance(1) == invalid_marker_variance);
}

TEST_CASE(
    "pooled Gaussian kernel matches the legacy transition under a fixed seed",
    "[bayes][kernel]")
{
    const auto model = make_model();
    const auto prior
        = gelex::compile(gelex::BayesRecipe<Method, mode_a>::defaults(), model);
    auto state = gelex::make_state(prior, model);
    auto kernel = gelex::make_kernel(prior);

    const auto legacy_model = make_model();
    const auto legacy_prior = make_legacy_prior(
        prior,
        {make_legacy_gaussian_prior(
            gelex::GeneticMode::A,
            prior.genetic().get<gelex::GeneticMode::A>().variance)});
    auto legacy_state = gelex::LegacyBayesState{legacy_model, legacy_prior};
    std::mt19937_64 rng{123};
    std::mt19937_64 legacy_rng{123};

    kernel.step(model, state, rng);
    auto legacy_kernel = gelex::Chain::make(
        legacy_model, legacy_prior, legacy_state, legacy_rng);
    legacy_kernel.step();

    const auto& genetic = state.genetic().get<gelex::GeneticMode::A>();
    const auto& legacy_block = std::get<gelex::bayes::SingleGeneticBlockState>(
        legacy_state.genetics().front());
    const auto& legacy_family
        = std::get<gelex::bayes::SingleSharedGaussianState>(
            legacy_block.prior_state());

    REQUIRE(state.fixed().coefficients.isApprox(legacy_state.fixed().coeffs));
    REQUIRE(genetic.coefficients.isApprox(legacy_block.state().coeffs));
    REQUIRE(state.residual().adjusted_response.isApprox(
        legacy_state.residual().y_adj));
    REQUIRE(genetic.family_state.variance == Approx(legacy_family.variance()));
    REQUIRE(
        state.residual().variance == Approx(legacy_state.residual().variance));
}

TEST_CASE(
    "independent AD pooled Gaussian kernel matches the legacy canonical mode "
    "order under a fixed seed",
    "[bayes][kernel]")
{
    const auto model = make_ad_model();
    const auto prior = gelex::compile(
        gelex::BayesRecipe<Method, mode_ad>::defaults(), model);
    auto state = gelex::make_state(prior, model);
    auto kernel = gelex::make_kernel(prior);

    const auto legacy_model = make_ad_model();
    const auto legacy_prior = make_legacy_prior(
        prior,
        {make_legacy_gaussian_prior(
             gelex::GeneticMode::A,
             prior.genetic().get<gelex::GeneticMode::A>().variance),
         make_legacy_gaussian_prior(
             gelex::GeneticMode::D,
             prior.genetic().get<gelex::GeneticMode::D>().variance)});
    auto legacy_state = gelex::LegacyBayesState{legacy_model, legacy_prior};
    std::mt19937_64 rng{123};
    std::mt19937_64 legacy_rng{123};

    kernel.step(model, state, rng);
    auto legacy_kernel = gelex::Chain::make(
        legacy_model, legacy_prior, legacy_state, legacy_rng);
    legacy_kernel.step();

    const auto& additive = state.genetic().get<gelex::GeneticMode::A>();
    const auto& dominance = state.genetic().get<gelex::GeneticMode::D>();
    const auto& legacy_additive
        = std::get<gelex::bayes::SingleGeneticBlockState>(
            legacy_state.genetics().at(0));
    const auto& legacy_dominance
        = std::get<gelex::bayes::SingleGeneticBlockState>(
            legacy_state.genetics().at(1));
    const auto& legacy_additive_family
        = std::get<gelex::bayes::SingleSharedGaussianState>(
            legacy_additive.prior_state());
    const auto& legacy_dominance_family
        = std::get<gelex::bayes::SingleSharedGaussianState>(
            legacy_dominance.prior_state());

    REQUIRE(state.fixed().coefficients.isApprox(legacy_state.fixed().coeffs));
    REQUIRE(additive.coefficients.isApprox(legacy_additive.state().coeffs));
    REQUIRE(dominance.coefficients.isApprox(legacy_dominance.state().coeffs));
    REQUIRE(additive.component_fitted_values.col(0).isApprox(
        legacy_additive.state().u));
    REQUIRE(dominance.component_fitted_values.col(0).isApprox(
        legacy_dominance.state().u));
    REQUIRE(state.residual().adjusted_response.isApprox(
        legacy_state.residual().y_adj));
    REQUIRE(
        additive.family_state.variance
        == Approx(legacy_additive_family.variance()));
    REQUIRE(
        dominance.family_state.variance
        == Approx(legacy_dominance_family.variance()));
    REQUIRE(
        state.residual().variance == Approx(legacy_state.residual().variance));
}

TEST_CASE(
    "unpooled Gaussian kernel matches the legacy transition under a fixed "
    "seed",
    "[bayes][kernel]")
{
    const auto model = make_model();
    const auto prior = gelex::compile(
        gelex::BayesRecipe<UnpooledMethod, mode_a>::defaults(), model);
    auto state = gelex::make_state(prior, model);
    auto kernel = gelex::make_kernel(prior);

    const auto legacy_model = make_model();
    const auto legacy_prior = make_legacy_prior(
        prior,
        {make_legacy_unpooled_gaussian_prior(
            gelex::GeneticMode::A,
            prior.genetic().get<gelex::GeneticMode::A>().variance)});
    auto legacy_state = gelex::LegacyBayesState{legacy_model, legacy_prior};
    std::mt19937_64 rng{123};
    std::mt19937_64 legacy_rng{123};

    kernel.step(model, state, rng);
    auto legacy_kernel = gelex::Chain::make(
        legacy_model, legacy_prior, legacy_state, legacy_rng);
    legacy_kernel.step();

    const auto& genetic = state.genetic().get<gelex::GeneticMode::A>();
    const auto& legacy_block = std::get<gelex::bayes::SingleGeneticBlockState>(
        legacy_state.genetics().front());
    const auto& legacy_family
        = std::get<gelex::bayes::SinglePerMarkerGaussianState>(
            legacy_block.prior_state());

    REQUIRE(state.fixed().coefficients.isApprox(legacy_state.fixed().coeffs));
    REQUIRE(genetic.coefficients.isApprox(legacy_block.state().coeffs));
    REQUIRE(genetic.component_fitted_values.col(0).isApprox(
        legacy_block.state().u));
    REQUIRE(state.residual().adjusted_response.isApprox(
        legacy_state.residual().y_adj));
    REQUIRE(genetic.family_state.variance.isApprox(legacy_family.variance()));
    REQUIRE(
        state.residual().variance == Approx(legacy_state.residual().variance));
}

TEST_CASE(
    "independent AD pooled Gaussian kernel rejects a missing mode before "
    "mutation",
    "[bayes][kernel]")
{
    const auto model = make_ad_model();
    const auto prior = gelex::compile(
        gelex::BayesRecipe<Method, mode_ad>::defaults(), model);
    auto state = gelex::make_state(prior, model);
    auto kernel = gelex::make_kernel(prior);
    const auto incomplete_model = make_model();
    std::mt19937_64 rng{123};
    const auto initial_rng = rng;

    REQUIRE_THROWS_AS(
        kernel.step(incomplete_model, state, rng), gelex::GelexException);

    REQUIRE(state.fixed().coefficients.isZero());
    REQUIRE(state.genetic().get<gelex::GeneticMode::A>().coefficients.isZero());
    REQUIRE(state.genetic().get<gelex::GeneticMode::D>().coefficients.isZero());
    REQUIRE(state.residual().adjusted_response.isApprox(model.phenotype()));
    REQUIRE(rng == initial_rng);
}

TEST_CASE(
    "pooled Gaussian kernel preserves random-effect and residual invariants",
    "[bayes][kernel]")
{
    const auto model = make_model_with_random();
    const auto prior = gelex::compile(
        gelex::BayesRecipe<Method, mode_a>{
            gelex::VarianceBudget{{.additive = 0.4, .random = 0.1}}},
        model);
    auto state = gelex::make_state(prior, model);
    auto kernel = gelex::make_kernel(prior);
    std::mt19937_64 rng{123};

    kernel.step(model, state, rng);

    const auto& random = state.random().front();
    const auto& genetic = state.genetic().get<gelex::GeneticMode::A>();
    const Eigen::VectorXd fixed_fitted
        = model.fixed().X * state.fixed().coefficients;
    const Eigen::VectorXd random_fitted
        = model.random().front().X * random.coefficients;
    const auto genetic_fitted = reconstruct_genetic_fitted(
        model.genetic(), gelex::GeneticMode::A, genetic.coefficients);

    REQUIRE(random.coefficients.allFinite());
    REQUIRE(random.variance > 0.0);
    REQUIRE((state.residual().adjusted_response + fixed_fitted + random_fitted
             + genetic_fitted)
                .isApprox(model.phenotype()));
}

TEST_CASE(
    "random-effect kernel matches the legacy transition under a fixed seed",
    "[bayes][kernel]")
{
    const auto model = make_model_with_random();
    const auto prior = gelex::compile(
        gelex::BayesRecipe<Method, mode_a>{
            gelex::VarianceBudget{{.additive = 0.4, .random = 0.1}}},
        model);
    auto state = gelex::make_state(prior, model);
    const auto& parameter = prior.random().front();
    gelex::bayes::RandomPrior legacy_prior{gelex::bayes::VarianceParameter{
        parameter.initial_value(),
        gelex::bayes::ScaledInvChiSqPrior{
            parameter.prior().degrees_of_freedom(),
            parameter.prior().scale()}}};
    std::vector<gelex::bayes::RandomState> legacy_states;
    legacy_states.emplace_back(model.random().front(), legacy_prior);
    gelex::bayes::ResidualState legacy_residual{
        .y_adj = state.residual().adjusted_response,
        .variance = state.residual().variance};
    std::mt19937_64 typed_rng{123};
    std::mt19937_64 legacy_rng{123};

    gelex::detail::RandomEffectKernel kernel{parameter};
    kernel.step(
        model.random().front(),
        state.random().front(),
        state.residual(),
        typed_rng);
    gelex::RandomStep legacy_kernel{
        legacy_prior,
        model.random(),
        legacy_states,
        legacy_residual,
        legacy_rng};
    legacy_kernel.step();

    REQUIRE(state.random().front().coefficients.isApprox(
        legacy_states.front().coeffs));
    REQUIRE(
        state.random().front().variance
        == Approx(legacy_states.front().variance));
    REQUIRE(state.residual().adjusted_response.isApprox(legacy_residual.y_adj));
    REQUIRE(typed_rng == legacy_rng);
}

TEST_CASE(
    "pooled spike-slab kernel preserves collapsed-state invariants",
    "[bayes][kernel]")
{
    const auto model = make_model();
    const auto prior = gelex::compile(
        gelex::BayesRecipe<PooledSpikeSlabMethod, mode_a>{
            gelex::IndependentTopology<mode_a, gelex::SpikeSlab>{
                gelex::SpikeSlab{.probability = 0.5}},
            gelex::VarianceBudget{{.additive = 0.4}}},
        model);
    auto state = gelex::make_state(prior, model);
    auto kernel = gelex::make_kernel(prior);
    std::mt19937_64 rng{123};

    kernel.step(model, state, rng);

    const auto& genetic = state.genetic().get<gelex::GeneticMode::A>();
    const auto& family = genetic.family_state;
    const auto reconstructed = reconstruct_genetic_fitted(
        model.genetic(), gelex::GeneticMode::A, genetic.coefficients);
    const Eigen::VectorXd fixed_fitted
        = model.fixed().X * state.fixed().coefficients;

    REQUIRE(genetic.component_fitted_values.col(0).isApprox(reconstructed));
    REQUIRE((state.residual().adjusted_response + fixed_fitted + reconstructed)
                .isApprox(model.phenotype()));
    REQUIRE(family.assignment.size() == model.genetic().cols());
    for (Eigen::Index marker = 0; marker < family.assignment.size(); ++marker)
    {
        REQUIRE(family.assignment(marker) <= 1);
        if (family.assignment(marker) == 0)
        {
            REQUIRE(genetic.coefficients(marker) == 0.0);
        }
    }
    REQUIRE(family.assignment(1) == 0);
    REQUIRE(genetic.coefficients(1) == 0.0);
    REQUIRE(family.variance > 0.0);
    REQUIRE(family.probability > 0.0);
    REQUIRE(family.probability < 1.0);
    REQUIRE(family.probability != 0.5);
}

TEST_CASE(
    "independent AD pooled spike-slab kernel preserves mode caches",
    "[bayes][kernel]")
{
    const auto model = make_ad_model();
    const auto prior = gelex::compile(
        gelex::BayesRecipe<PooledSpikeSlabMethod, mode_ad>{
            gelex::IndependentTopology<mode_ad, gelex::SpikeSlab>{
                gelex::SpikeSlab{.probability = 0.5},
                gelex::SpikeSlab{.probability = 0.5}},
            gelex::VarianceBudget{{.additive = 0.4, .dominance = 0.1}}},
        model);
    auto state = gelex::make_state(prior, model);
    auto kernel = gelex::make_kernel(prior);
    std::mt19937_64 rng{123};

    kernel.step(model, state, rng);

    const auto& additive = state.genetic().get<gelex::GeneticMode::A>();
    const auto& dominance = state.genetic().get<gelex::GeneticMode::D>();
    const auto additive_fitted = reconstruct_genetic_fitted(
        model.genetic(), gelex::GeneticMode::A, additive.coefficients);
    const auto dominance_fitted = reconstruct_genetic_fitted(
        model.genetic(), gelex::GeneticMode::D, dominance.coefficients);
    const Eigen::VectorXd fixed_fitted
        = model.fixed().X * state.fixed().coefficients;

    REQUIRE(additive.component_fitted_values.col(0).isApprox(additive_fitted));
    REQUIRE(
        dominance.component_fitted_values.col(0).isApprox(dominance_fitted));
    REQUIRE((state.residual().adjusted_response + fixed_fitted + additive_fitted
             + dominance_fitted)
                .isApprox(model.phenotype()));
    for (Eigen::Index marker = 0; marker < model.genetic().cols(); ++marker)
    {
        if (additive.family_state.assignment(marker) == 0)
        {
            REQUIRE(additive.coefficients(marker) == 0.0);
        }
        if (dominance.family_state.assignment(marker) == 0)
        {
            REQUIRE(dominance.coefficients(marker) == 0.0);
        }
    }
}

TEST_CASE(
    "unpooled spike-slab kernel samples inactive variances and keeps a fixed "
    "probability",
    "[bayes][kernel]")
{
    const auto model = make_model();
    constexpr double fixed_probability = 1.0e-12;
    const auto prior = gelex::compile(
        gelex::BayesRecipe<FixedUnpooledSpikeSlabMethod, mode_a>{
            gelex::IndependentTopology<mode_a, gelex::SpikeSlab>{
                gelex::SpikeSlab{.probability = fixed_probability}},
            gelex::VarianceBudget{{.additive = 0.4}}},
        model);
    auto state = gelex::make_state(prior, model);
    auto kernel = gelex::make_kernel(prior);
    auto& initial_family
        = state.genetic().get<gelex::GeneticMode::A>().family_state;
    const Eigen::VectorXd initial_variance = initial_family.variance;
    std::mt19937_64 rng{123};

    kernel.step(model, state, rng);

    const auto& genetic = state.genetic().get<gelex::GeneticMode::A>();
    const auto& family = genetic.family_state;
    const Eigen::VectorXd fixed_fitted
        = model.fixed().X * state.fixed().coefficients;

    REQUIRE(genetic.coefficients.isZero());
    REQUIRE(family.assignment.isZero());
    REQUIRE(genetic.component_fitted_values.col(0).isZero());
    REQUIRE((state.residual().adjusted_response + fixed_fitted)
                .isApprox(model.phenotype()));
    REQUIRE((family.variance.array() > 0.0).all());
    REQUIRE(family.variance(0) != initial_variance(0));
    REQUIRE(family.variance(1) == initial_variance(1));
    REQUIRE(family.probability == fixed_probability);
}

TEST_CASE(
    "scaled-mixture kernel samples classes and probabilities",
    "[bayes][kernel]")
{
    const auto model = make_model();
    constexpr std::array probabilities{0.05, 0.1, 0.15, 0.2, 0.5};
    const auto prior = gelex::compile(
        gelex::BayesRecipe<SampledScaledMixtureMethod, mode_a>{
            gelex::IndependentTopology<mode_a, gelex::ScaledMixture>{
                gelex::ScaledMixture{.probabilities = probabilities}},
            gelex::VarianceBudget{{.additive = 0.4}}},
        model);
    auto state = gelex::make_state(prior, model);
    auto kernel = gelex::make_kernel(prior);
    std::mt19937_64 rng{123};

    kernel.step(model, state, rng);

    const auto& genetic = state.genetic().get<gelex::GeneticMode::A>();
    const auto& family = genetic.family_state;
    const auto fitted = reconstruct_genetic_fitted(
        model.genetic(), gelex::GeneticMode::A, genetic.coefficients);
    const auto component_fitted = reconstruct_component_fitted<
        gelex::ScaledMixturePrior<>::component_layout>(
        model.genetic(),
        gelex::GeneticMode::A,
        genetic.coefficients,
        family.assignment);
    const Eigen::VectorXd fixed_fitted
        = model.fixed().X * state.fixed().coefficients;
    REQUIRE((state.residual().adjusted_response + fixed_fitted + fitted)
                .isApprox(model.phenotype()));
    REQUIRE(genetic.component_fitted_values.isApprox(component_fitted));
    REQUIRE(family.assignment(1) == 0);
    REQUIRE(genetic.coefficients(1) == 0.0);
    for (Eigen::Index marker = 0; marker < family.assignment.size(); ++marker)
    {
        REQUIRE(
            family.assignment(marker)
            < gelex::ScaledMixturePrior<>::class_count);
        if (family.assignment(marker) == 0)
        {
            REQUIRE(genetic.coefficients(marker) == 0.0);
        }
    }
    double probability_sum = 0.0;
    for (const double probability : family.probabilities)
    {
        REQUIRE(probability > 0.0);
        probability_sum += probability;
    }
    REQUIRE(probability_sum == Approx(1.0));
    REQUIRE(family.probabilities != probabilities);
    REQUIRE(family.variance > 0.0);
}

TEST_CASE(
    "fixed scaled-mixture kernel omits probability updates",
    "[bayes][kernel]")
{
    const auto model = make_model();
    constexpr std::array probabilities{
        1.0e-12, 1.0e-12, 1.0e-12, 1.0e-12, 1.0 - 4.0e-12};
    const auto prior = gelex::compile(
        gelex::BayesRecipe<FixedScaledMixtureMethod, mode_a>{
            gelex::IndependentTopology<mode_a, gelex::ScaledMixture>{
                gelex::ScaledMixture{.probabilities = probabilities}},
            gelex::VarianceBudget{{.additive = 0.4}}},
        model);
    auto state = gelex::make_state(prior, model);
    auto kernel = gelex::make_kernel(prior);
    std::mt19937_64 rng{123};

    kernel.step(model, state, rng);

    const auto& genetic = state.genetic().get<gelex::GeneticMode::A>();
    const auto& family = genetic.family_state;
    const auto reconstructed = reconstruct_genetic_fitted(
        model.genetic(), gelex::GeneticMode::A, genetic.coefficients);
    const auto component_fitted = reconstruct_component_fitted<
        gelex::ScaledMixturePrior<>::component_layout>(
        model.genetic(),
        gelex::GeneticMode::A,
        genetic.coefficients,
        family.assignment);
    const Eigen::VectorXd fixed_fitted
        = model.fixed().X * state.fixed().coefficients;

    REQUIRE(family.probabilities == probabilities);
    REQUIRE(genetic.component_fitted_values.isApprox(component_fitted));
    REQUIRE((state.residual().adjusted_response + fixed_fitted + reconstructed)
                .isApprox(model.phenotype()));
    REQUIRE(family.variance > 0.0);
}

TEST_CASE(
    "joint spike-slab kernel preserves blocked-state invariants",
    "[bayes][kernel]")
{
    const auto model = make_joint_model();
    constexpr std::array probabilities{0.25, 0.25, 0.25, 0.25};
    constexpr double positive_probability = 0.6;
    const auto prior = gelex::compile(
        gelex::BayesRecipe<SampledJointSpikeSlabMethod, mode_ad>{
            gelex::JointSpikeSlab{
                .probabilities = probabilities,
                .positive_probability = positive_probability},
            gelex::VarianceBudget{{.additive = 0.4, .dominance = 0.1}}},
        model);
    auto state = gelex::make_state(prior, model);
    auto kernel = gelex::make_kernel(prior);
    std::mt19937_64 rng{123};

    REQUIRE(
        model.genetic()
            .projection(gelex::GeneticMode::A)
            .col_covariance(model.genetic().projection(gelex::GeneticMode::D))
            .isZero(1.0e-12));

    kernel.step(model, state, rng);

    const auto& mode_states = state.genetic().mode_values();
    const auto& additive = mode_states.get<gelex::GeneticMode::A>();
    const auto& dominance = mode_states.get<gelex::GeneticMode::D>();
    const auto& joint = state.genetic().joint();
    const auto additive_fitted = reconstruct_genetic_fitted(
        model.genetic(), gelex::GeneticMode::A, additive.coefficients);
    const auto dominance_fitted = reconstruct_genetic_fitted(
        model.genetic(), gelex::GeneticMode::D, dominance.coefficients);
    const auto additive_component_fitted = reconstruct_component_fitted<
        gelex::JointSpikeSlabPrior<>::component_layout>(
        model.genetic(),
        gelex::GeneticMode::A,
        additive.coefficients,
        joint.assignment);
    const auto dominance_component_fitted = reconstruct_component_fitted<
        gelex::JointSpikeSlabPrior<>::component_layout>(
        model.genetic(),
        gelex::GeneticMode::D,
        dominance.coefficients,
        joint.assignment);
    const Eigen::VectorXd fixed_fitted
        = model.fixed().X * state.fixed().coefficients;

    REQUIRE((state.residual().adjusted_response + fixed_fitted + additive_fitted
             + dominance_fitted)
                .isApprox(model.phenotype()));
    REQUIRE(
        additive.component_fitted_values.isApprox(additive_component_fitted));
    REQUIRE(
        dominance.component_fitted_values.isApprox(dominance_component_fitted));
    for (Eigen::Index marker = 0; marker < joint.assignment.size(); ++marker)
    {
        const auto class_index = joint.assignment(marker);
        REQUIRE(class_index < gelex::JointSpikeSlabPrior<>::class_count);
        const bool additive_active = class_index == 1 || class_index == 3;
        const bool dominance_active = class_index == 2 || class_index == 3;
        if (!additive_active)
        {
            REQUIRE(additive.coefficients(marker) == 0.0);
        }
        if (dominance_active)
        {
            REQUIRE(
                (dominance.coefficients(marker) > 0.0)
                == (joint.dominance_sign(marker) == 1));
        }
        else
        {
            REQUIRE(dominance.coefficients(marker) == 0.0);
            REQUIRE(joint.dominance_sign(marker) == 0);
        }
    }
    REQUIRE(joint.assignment(1) == 0);
    REQUIRE(additive.coefficients(1) == 0.0);
    REQUIRE(dominance.coefficients(1) == 0.0);
    REQUIRE(additive.family_state.variance > 0.0);
    REQUIRE(dominance.family_state.variance > 0.0);

    double probability_sum = 0.0;
    for (const double probability : joint.probabilities)
    {
        REQUIRE(probability > 0.0);
        probability_sum += probability;
    }
    REQUIRE(probability_sum == Approx(1.0));
    REQUIRE(joint.probabilities != probabilities);
    REQUIRE(joint.positive_probability > 0.0);
    REQUIRE(joint.positive_probability < 1.0);
    REQUIRE(joint.positive_probability != positive_probability);
}

TEST_CASE(
    "fixed joint spike-slab kernel omits probability updates",
    "[bayes][kernel]")
{
    const auto model = make_joint_model();
    constexpr std::array probabilities{0.1, 0.2, 0.3, 0.4};
    constexpr double positive_probability = 0.75;
    const auto prior = gelex::compile(
        gelex::BayesRecipe<FixedJointSpikeSlabMethod, mode_ad>{
            gelex::JointSpikeSlab{
                .probabilities = probabilities,
                .positive_probability = positive_probability},
            gelex::VarianceBudget{{.additive = 0.4, .dominance = 0.1}}},
        model);
    auto state = gelex::make_state(prior, model);
    auto kernel = gelex::make_kernel(prior);
    std::mt19937_64 rng{123};

    kernel.step(model, state, rng);

    const auto& mode_states = state.genetic().mode_values();
    const auto& additive = mode_states.get<gelex::GeneticMode::A>();
    const auto& dominance = mode_states.get<gelex::GeneticMode::D>();
    const auto& joint = state.genetic().joint();
    const auto additive_fitted = reconstruct_genetic_fitted(
        model.genetic(), gelex::GeneticMode::A, additive.coefficients);
    const auto dominance_fitted = reconstruct_genetic_fitted(
        model.genetic(), gelex::GeneticMode::D, dominance.coefficients);
    const Eigen::VectorXd fixed_fitted
        = model.fixed().X * state.fixed().coefficients;

    REQUIRE(joint.probabilities == probabilities);
    REQUIRE(joint.positive_probability == positive_probability);
    REQUIRE((state.residual().adjusted_response + fixed_fitted + additive_fitted
             + dominance_fitted)
                .isApprox(model.phenotype()));
    REQUIRE(additive.family_state.variance > 0.0);
    REQUIRE(dominance.family_state.variance > 0.0);
}

TEST_CASE(
    "joint spike-slab kernel matches the NOIA legacy transition under a fixed "
    "seed",
    "[bayes][kernel]")
{
    constexpr std::array probabilities{0.25, 0.25, 0.25, 0.25};
    constexpr double positive_probability = 0.6;
    const auto model = make_joint_model();
    const auto prior = gelex::compile(
        gelex::BayesRecipe<SampledJointSpikeSlabMethod, mode_ad>{
            gelex::JointSpikeSlab{
                .probabilities = probabilities,
                .positive_probability = positive_probability},
            gelex::VarianceBudget{{.additive = 0.4, .dominance = 0.1}}},
        model);
    auto state = gelex::make_state(prior, model);
    auto kernel = gelex::make_kernel(prior);

    const auto legacy_model = make_joint_model();
    const auto legacy_prior
        = make_legacy_prior(prior, {make_legacy_joint_prior(prior.genetic())});
    auto legacy_state = gelex::LegacyBayesState{legacy_model, legacy_prior};
    std::mt19937_64 rng{123};
    std::mt19937_64 legacy_rng{123};

    kernel.step(model, state, rng);
    auto legacy_kernel = gelex::Chain::make(
        legacy_model, legacy_prior, legacy_state, legacy_rng);
    legacy_kernel.step();

    const auto& mode_states = state.genetic().mode_values();
    const auto& additive = mode_states.get<gelex::GeneticMode::A>();
    const auto& dominance = mode_states.get<gelex::GeneticMode::D>();
    const auto& joint = state.genetic().joint();
    const auto& legacy_block = std::get<gelex::bayes::JointGeneticBlockState>(
        legacy_state.genetics().front());
    const auto& legacy_joint
        = std::get<gelex::bayes::JointHalfNormalMixtureState>(
            legacy_block.prior_state());

    REQUIRE(state.fixed().coefficients.isApprox(legacy_state.fixed().coeffs));
    REQUIRE(additive.coefficients.isApprox(
        legacy_block.state(gelex::GeneticMode::A).coeffs));
    REQUIRE(dominance.coefficients.isApprox(
        legacy_block.state(gelex::GeneticMode::D).coeffs));
    REQUIRE(state.residual().adjusted_response.isApprox(
        legacy_state.residual().y_adj));
    REQUIRE(
        additive.family_state.variance
        == Approx(legacy_joint.variance(gelex::GeneticMode::A)));
    REQUIRE(
        dominance.family_state.variance
        == Approx(legacy_joint.variance(gelex::GeneticMode::D)));
    REQUIRE(
        state.residual().variance == Approx(legacy_state.residual().variance));
    REQUIRE(
        joint.positive_probability
        == Approx(legacy_joint.dominance_sign().positive_probability));
    for (Eigen::Index marker = 0; marker < joint.assignment.size(); ++marker)
    {
        REQUIRE(joint.assignment(marker) == legacy_joint.assignment()(marker));
        REQUIRE(
            joint.dominance_sign(marker)
            == legacy_joint.dominance_sign().sign(marker));
    }
    for (std::size_t class_index = 0; class_index < probabilities.size();
         ++class_index)
    {
        REQUIRE(
            joint.probabilities.at(class_index)
            == Approx(legacy_joint.proportion()(
                static_cast<Eigen::Index>(class_index))));
    }
    REQUIRE(rng == legacy_rng);
}
