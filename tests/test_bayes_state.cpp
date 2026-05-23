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
#include <memory>
#include <span>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

#include <Eigen/Core>
#include <catch2/catch_test_macros.hpp>

#include "gelex/data/genotype/genotype.h"
#include "gelex/exception.h"
#include "gelex/infra/record_visitor.h"
#include "gelex/model/bayes/gaussian_prior.h"
#include "gelex/model/bayes/genetic_prior.h"
#include "gelex/model/bayes/model.h"
#include "gelex/model/bayes/prior.h"
#include "gelex/model/bayes/prior_specs.h"
#include "gelex/model/bayes/state.h"
#include "gelex/model/bayes/state_capabilities.h"
#include "gelex/model/bayes/recipe.h"
#include "gelex/model/bayes/recipe_options.h"
#include "gelex/types/fixed_effects.h"
#include "gelex/types/genetic_effect_type.h"
#include "genotype_fixture.h"

using gelex::BayesModel;
using gelex::BayesState;
using gelex::FixedEffect;
using gelex::GeneticMode;
using gelex::GelexException;
using gelex::bayes::BayesPrior;
using gelex::bayes::BayesRecipe;
using gelex::bayes::BayesRecipeConfig;
using gelex::bayes::BayesRecipePreset;
using gelex::bayes::ComponentStateCap;
using gelex::bayes::GaussianPrior;
using gelex::bayes::GeneticPrior;
using gelex::bayes::JointMixtureGaussianPrior;
using gelex::bayes::MarkerVarianceScope;
using gelex::bayes::MarkerVarianceSpec;
using gelex::bayes::ProportionSpec;
using gelex::bayes::ProportionStateCap;
using gelex::bayes::ProportionUpdate;
using gelex::bayes::ScaledInvChiSqPrior;
using gelex::bayes::ScaledMixtureGaussianPrior;
using gelex::bayes::SpikeSlabGaussianPrior;
using gelex::bayes::VarianceSpec;
using gelex::bayes::VarianceStateCap;

namespace
{

constexpr Eigen::Index kNumIndividuals = 4;
constexpr Eigen::Index kNumMarkers = 3;

auto make_genotype(Eigen::MatrixXd data) -> gelex::genotype::Genotype
{
    const Eigen::Index cols = data.cols();
    auto mean = data.colwise().mean().transpose().eval();
    auto stddev = Eigen::VectorXd::Ones(cols);
    return gelex::test::GenotypeBuilder::build(
        std::move(data), std::move(mean), std::move(stddev));
}

auto make_genetic_data(GeneticMode mode) -> Eigen::MatrixXd
{
    Eigen::MatrixXd data{
        {0.0, 1.0, 2.0},
        {1.0, 2.0, 0.0},
        {2.0, 0.0, 1.0},
        {1.0, 1.0, 2.0},
    };
    if (mode == GeneticMode::D)
    {
        data.array() += 0.5;
    }
    return data;
}

auto make_model(
    std::span<const GeneticMode> modes,
    bool include_random = false) -> BayesModel
{
    Eigen::VectorXd phenotype{{1.0, 2.0, 4.0, 8.0}};
    auto fixed = FixedEffect::build(kNumIndividuals);

    std::vector<gelex::bayes::RandomEffect> random;
    if (include_random)
    {
        random.emplace_back(
            "batch",
            std::vector<std::string>{"a", "b"},
            Eigen::MatrixXd{{1.0, 0.0}, {1.0, 0.0}, {0.0, 1.0}, {0.0, 1.0}});
    }

    std::vector<gelex::bayes::GeneticEffect> genetics;
    for (const auto mode : modes)
    {
        genetics.emplace_back(mode, make_genotype(make_genetic_data(mode)));
    }

    return BayesModel(
        std::move(phenotype),
        std::move(fixed),
        std::move(random),
        std::move(genetics));
}

auto make_variance(double initial_value = 1.0) -> VarianceSpec
{
    return VarianceSpec(initial_value, ScaledInvChiSqPrior{4.0, 1.0});
}

auto make_marker_variance(
    MarkerVarianceScope scope = MarkerVarianceScope::per_marker,
    double initial_value = 1.0) -> MarkerVarianceSpec
{
    return MarkerVarianceSpec{scope, make_variance(initial_value)};
}

auto make_proportion_2() -> ProportionSpec
{
    return ProportionSpec{
        Eigen::VectorXd{{0.9, 0.1}},
        Eigen::VectorXd{{1.0, 1.0}},
        ProportionUpdate::fixed,
    };
}

auto make_proportion_4() -> ProportionSpec
{
    return ProportionSpec{
        Eigen::VectorXd{{0.91, 0.03, 0.03, 0.03}},
        Eigen::VectorXd{{1.0, 1.0, 1.0, 1.0}},
        ProportionUpdate::sampled,
    };
}

auto make_prior(
    std::vector<std::unique_ptr<GeneticPrior>> genetics) -> BayesPrior
{
    return BayesPrior(make_variance(0.25), std::move(genetics), make_variance(1.5));
}

class PathCollector final : public gelex::infra::RecordSink
{
   public:
    auto visit(std::string_view path, const Eigen::Ref<const Eigen::VectorXf>&)
        -> void override
    {
        paths.emplace_back(path);
    }

    auto visit(std::string_view path, const Eigen::Ref<const Eigen::VectorXd>&)
        -> void override
    {
        paths.emplace_back(path);
    }

    auto visit(std::string_view path, const Eigen::Ref<const Eigen::VectorXi>&)
        -> void override
    {
        paths.emplace_back(path);
    }

    auto visit(std::string_view path, const double&) -> void override
    {
        paths.emplace_back(path);
    }

    auto count(std::string_view path) const -> std::size_t
    {
        return static_cast<std::size_t>(
            std::ranges::count(paths, std::string{path}));
    }

    std::vector<std::string> paths;
};

}  // namespace

TEST_CASE("BayesRecipe prior constructs BayesState", "[bayes_state]")
{
    constexpr std::array modes{GeneticMode::A};
    auto model = make_model(modes, true);

    BayesRecipeConfig config;
    config.modes = {GeneticMode::A};
    BayesRecipe recipe(BayesRecipePreset::RR, config);
    auto prior = recipe.make_prior(model);

    BayesState state(model, prior);

    REQUIRE(state.fixed().coeffs.size() == 1);
    REQUIRE(state.fixed().coeffs.isZero());
    REQUIRE(state.random().size() == 1);
    REQUIRE(state.random()[0].coeffs.size() == 2);
    REQUIRE(state.genetics().size() == 1);
    REQUIRE(state.genetic(GeneticMode::A) != nullptr);
    REQUIRE(state.residual().y_adj.isApprox(model.phenotype()));
    REQUIRE(state.residual().variance > 0.0);
}

TEST_CASE("BayesState initializes fixed random and residual state", "[bayes_state]")
{
    constexpr std::array modes{GeneticMode::A};
    auto model = make_model(modes, true);

    std::vector<std::unique_ptr<GeneticPrior>> genetics;
    genetics.push_back(std::make_unique<GaussianPrior>(
        GeneticMode::A, make_marker_variance()));
    auto prior = make_prior(std::move(genetics));

    BayesState state(model, prior);

    REQUIRE(state.fixed().coeffs.isApprox(Eigen::VectorXd::Zero(1)));
    REQUIRE(state.random().size() == 1);
    REQUIRE(state.random()[0].coeffs.isApprox(Eigen::VectorXd::Zero(2)));
    REQUIRE(state.random()[0].variance == 0.25);
    REQUIRE(state.residual().y_adj.isApprox(model.phenotype()));
    REQUIRE(state.residual().variance == 1.5);
}

TEST_CASE("BayesState creates genetic prior states by capability", "[bayes_state]")
{
    constexpr std::array modes{GeneticMode::A, GeneticMode::D};
    auto model = make_model(modes);

    SECTION("Gaussian")
    {
        std::vector<std::unique_ptr<GeneticPrior>> genetics;
        genetics.push_back(std::make_unique<GaussianPrior>(
            GeneticMode::A, make_marker_variance()));
        auto prior = make_prior(std::move(genetics));

        BayesState state(model, prior);
        auto* genetic = state.genetic(GeneticMode::A);
        auto* block = state.genetic_block_for(GeneticMode::A);

        REQUIRE(genetic != nullptr);
        REQUIRE(block != nullptr);
        REQUIRE(genetic->coeffs.size() == kNumMarkers);
        REQUIRE(genetic->u.size() == kNumIndividuals);
        REQUIRE(state.genetic_blocks().size() == 1);
        REQUIRE(block->slot(GeneticMode::A) == 0);
        REQUIRE(block->prior_state().query<VarianceStateCap>() != nullptr);
        REQUIRE(block->prior_state().query<ProportionStateCap>() == nullptr);
    }

    SECTION("SpikeSlab")
    {
        std::vector<std::unique_ptr<GeneticPrior>> genetics;
        genetics.push_back(std::make_unique<SpikeSlabGaussianPrior>(
            GeneticMode::A, make_marker_variance(), make_proportion_2()));
        auto prior = make_prior(std::move(genetics));

        BayesState state(model, prior);
        auto* genetic = state.genetic(GeneticMode::A);
        auto* block = state.genetic_block_for(GeneticMode::A);
        REQUIRE(genetic != nullptr);
        REQUIRE(block != nullptr);

        auto& proportion = block->prior_state().require<ProportionStateCap>();

        REQUIRE(proportion.proportion().size() == 1);
        REQUIRE(proportion.proportion()[0].assignment.size() == kNumMarkers);
        REQUIRE(proportion.proportion()[0].count.isApprox(Eigen::VectorXi{{3, 0}}));
    }

    SECTION("ScaledMixture")
    {
        std::vector<std::unique_ptr<GeneticPrior>> genetics;
        genetics.push_back(std::make_unique<ScaledMixtureGaussianPrior>(
            GeneticMode::D,
            make_marker_variance(MarkerVarianceScope::per_effect),
            Eigen::VectorXd{{0.0, 0.1, 1.0}},
            ProportionSpec{
                Eigen::VectorXd{{0.8, 0.1, 0.1}},
                Eigen::VectorXd{{1.0, 1.0, 1.0}},
                ProportionUpdate::fixed}));
        auto prior = make_prior(std::move(genetics));

        BayesState state(model, prior);
        auto* genetic = state.genetic(GeneticMode::D);
        auto* block = state.genetic_block_for(GeneticMode::D);

        REQUIRE(genetic != nullptr);
        REQUIRE(block != nullptr);
        REQUIRE(block->prior_state().query<VarianceStateCap>() != nullptr);
        REQUIRE(block->prior_state().query<ProportionStateCap>() != nullptr);
        REQUIRE(block->prior_state().query<ComponentStateCap>() != nullptr);
    }
}

TEST_CASE("BayesState shares joint prior state across modes", "[bayes_state]")
{
    constexpr std::array modes{GeneticMode::A, GeneticMode::D};
    auto model = make_model(modes);

    std::vector<std::unique_ptr<GeneticPrior>> genetics;
    genetics.push_back(std::make_unique<JointMixtureGaussianPrior>(
        std::array{GeneticMode::A, GeneticMode::D},
        std::array{
            make_marker_variance(MarkerVarianceScope::per_effect),
            make_marker_variance(MarkerVarianceScope::per_effect)},
        make_proportion_4()));
    auto prior = make_prior(std::move(genetics));

    BayesState state(model, prior);
    auto* additive = state.genetic(GeneticMode::A);
    auto* dominance = state.genetic(GeneticMode::D);
    auto* additive_block = state.genetic_block_for(GeneticMode::A);
    auto* dominance_block = state.genetic_block_for(GeneticMode::D);

    REQUIRE(additive != nullptr);
    REQUIRE(dominance != nullptr);
    REQUIRE(additive_block != nullptr);
    REQUIRE(dominance_block != nullptr);
    REQUIRE(additive_block == dominance_block);
    REQUIRE(state.genetic_blocks().size() == 1);
    REQUIRE(additive_block->modes().size() == 2);
    REQUIRE(additive_block->genetic_indices().size() == 2);
    REQUIRE(additive_block->slot(GeneticMode::A) == 0);
    REQUIRE(additive_block->slot(GeneticMode::D) == 1);
    REQUIRE(&state.genetics()[additive_block->genetic_indices()[0]] == additive);
    REQUIRE(&state.genetics()[additive_block->genetic_indices()[1]] == dominance);
    REQUIRE(additive->coeffs.data() != dominance->coeffs.data());
    REQUIRE(additive->u.data() != dominance->u.data());

    auto& variance = additive_block->prior_state().require<VarianceStateCap>();
    auto& proportion
        = additive_block->prior_state().require<ProportionStateCap>();
    REQUIRE(variance.variance().size() == 2);
    REQUIRE(proportion.proportion().size() == 1);
    REQUIRE(proportion.proportion()[0].value.size() == 4);
}

TEST_CASE("BayesState rejects prior mode missing from model", "[bayes_state]")
{
    constexpr std::array modes{GeneticMode::A};
    auto model = make_model(modes);

    std::vector<std::unique_ptr<GeneticPrior>> genetics;
    genetics.push_back(std::make_unique<GaussianPrior>(
        GeneticMode::D, make_marker_variance()));
    auto prior = make_prior(std::move(genetics));

    REQUIRE_THROWS_AS(BayesState(model, prior), GelexException);
}

TEST_CASE("BayesState computes heritability from runtime variance", "[bayes_state]")
{
    constexpr std::array modes{GeneticMode::A, GeneticMode::D};
    auto model = make_model(modes, true);

    std::vector<std::unique_ptr<GeneticPrior>> genetics;
    genetics.push_back(std::make_unique<GaussianPrior>(
        GeneticMode::A, make_marker_variance()));
    genetics.push_back(std::make_unique<GaussianPrior>(
        GeneticMode::D, make_marker_variance()));
    auto prior = make_prior(std::move(genetics));

    BayesState state(model, prior);
    state.random()[0].variance = 2.0;
    state.genetic(GeneticMode::A)->variance = 3.0;
    state.genetic(GeneticMode::D)->variance = 5.0;
    state.residual().variance = 10.0;

    state.compute_heritability();

    REQUIRE(std::abs(state.genetic(GeneticMode::A)->heritability - 0.15) < 1e-12);
    REQUIRE(std::abs(state.genetic(GeneticMode::D)->heritability - 0.25) < 1e-12);
}

TEST_CASE("BayesState visits sample records through state schema", "[bayes_state]")
{
    constexpr std::array modes{GeneticMode::A};
    auto model = make_model(modes, true);

    std::vector<std::unique_ptr<GeneticPrior>> genetics;
    genetics.push_back(std::make_unique<SpikeSlabGaussianPrior>(
        GeneticMode::A, make_marker_variance(), make_proportion_2()));
    auto prior = make_prior(std::move(genetics));

    BayesState state(model, prior);
    state.genetic(GeneticMode::A)->variance = 2.0;
    state.compute_heritability();

    PathCollector sink;
    state.visit_records(gelex::bayes::StateRecordSet::sample, sink);

    REQUIRE(sink.count("fixed/0/coeffs") == 1);
    REQUIRE(sink.count("random/0/coeffs") == 1);
    REQUIRE(sink.count("random/0/variance") == 1);
    REQUIRE(sink.count("genetic/0/coeffs") == 1);
    REQUIRE(sink.count("genetic/0/variance") == 1);
    REQUIRE(sink.count("genetic/0/heritability") == 1);
    REQUIRE(sink.count("genetic_block/0/prior_state/variance/0/value") == 1);
    REQUIRE(sink.count("genetic_block/0/prior_state/proportion/0/assignment") == 1);
    REQUIRE(sink.count("residual/0/variance") == 1);
}

TEST_CASE("BayesState visits joint prior state once for checkpoint", "[bayes_state]")
{
    constexpr std::array modes{GeneticMode::A, GeneticMode::D};
    auto model = make_model(modes);

    std::vector<std::unique_ptr<GeneticPrior>> genetics;
    genetics.push_back(std::make_unique<JointMixtureGaussianPrior>(
        std::array{GeneticMode::A, GeneticMode::D},
        std::array{
            make_marker_variance(MarkerVarianceScope::per_effect),
            make_marker_variance(MarkerVarianceScope::per_effect)},
        make_proportion_4()));
    auto prior = make_prior(std::move(genetics));

    BayesState state(model, prior);
    PathCollector sink;
    state.visit_records(gelex::bayes::StateRecordSet::checkpoint, sink);

    REQUIRE(sink.count("genetic/0/coeffs") == 1);
    REQUIRE(sink.count("genetic/1/coeffs") == 1);
    REQUIRE(sink.count("genetic/0/u") == 1);
    REQUIRE(sink.count("genetic/1/u") == 1);
    REQUIRE(sink.count("genetic_block/0/prior_state/variance/0/value") == 1);
    REQUIRE(sink.count("genetic_block/0/prior_state/variance/1/value") == 1);
    REQUIRE(sink.count("genetic_block/0/prior_state/proportion/0/assignment") == 1);
    REQUIRE(sink.count("genetic_block/0/prior_state/proportion/0/count") == 1);
    REQUIRE(sink.count("genetic_block/0/prior_state/proportion/0/value") == 1);
    REQUIRE(sink.count("genetic_block/1/prior_state/proportion/0/value") == 0);
}
