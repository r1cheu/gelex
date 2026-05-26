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
#include <ranges>
#include <span>
#include <string>
#include <string_view>
#include <utility>
#include <variant>
#include <vector>

#include <Eigen/Core>
#include <catch2/catch_test_macros.hpp>

#include "gelex/data/genotype/genotype.h"
#include "gelex/exception.h"
#include "gelex/infra/field_visitor.h"
#include "gelex/infra/record_visitor.h"
#include "gelex/model/bayes/gaussian_prior.h"
#include "gelex/model/bayes/genetic_prior.h"
#include "gelex/model/bayes/genetic_priors/gaussian.h"
#include "gelex/model/bayes/model.h"
#include "gelex/model/bayes/prior.h"
#include "gelex/model/bayes/prior_parameters.h"
#include "gelex/model/bayes/recipe.h"
#include "gelex/model/bayes/recipe_options.h"
#include "gelex/model/bayes/state.h"
#include "gelex/model/bayes/state_capabilities.h"
#include "gelex/types/fixed_effects.h"
#include "gelex/types/genetic_effect_type.h"
#include "genotype_fixture.h"

using gelex::BayesModel;
using gelex::BayesState;
using gelex::BayesStateV2;
using gelex::FixedEffect;
using gelex::GelexException;
using gelex::GeneticMode;
using gelex::bayes::BayesPrior;
using gelex::bayes::BayesPriorV2;
using gelex::bayes::BayesRecipe;
using gelex::bayes::BayesRecipeConfig;
using gelex::bayes::BayesRecipePreset;
using gelex::bayes::ComponentStateCap;
using gelex::bayes::DirichletPrior;
using gelex::bayes::GaussianPrior;
using gelex::bayes::GeneticPriorBlockV2;
using gelex::bayes::GeneticPriorBlockState;
using gelex::bayes::GeneticPrior;
using gelex::bayes::JointGaussianMixturePrior;
using gelex::bayes::JointGeneticBlockState;
using gelex::bayes::JointMixtureGaussianPrior;
using gelex::bayes::JointVarianceStateCap;
using gelex::bayes::MarkerVariance;
using gelex::bayes::MarkerVarianceLayout;
using gelex::bayes::MixtureProportion;
using gelex::bayes::ProportionStateCap;
using gelex::bayes::RandomPrior;
using gelex::bayes::ResidualPrior;
using gelex::bayes::ScaledInvChiSqPrior;
using gelex::bayes::ScaledMixtureGaussianPrior;
using gelex::bayes::SimplexParameter;
using gelex::bayes::SingleGaussianPrior;
using gelex::bayes::SingleGeneticBlockState;
using gelex::bayes::SingleVarianceStateCap;
using gelex::bayes::SpikeSlabGaussianPrior;
using gelex::bayes::UpdatePolicy;
using gelex::bayes::VarianceParameter;
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

auto make_model(std::span<const GeneticMode> modes, bool include_random = false)
    -> BayesModel
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

auto make_variance(double initial_value = 1.0) -> VarianceParameter
{
    return VarianceParameter(initial_value, ScaledInvChiSqPrior{4.0, 1.0});
}

auto make_random_prior(double initial_value = 1.0) -> RandomPrior
{
    return RandomPrior{make_variance(initial_value)};
}

auto make_residual_prior(double initial_value = 1.0) -> ResidualPrior
{
    return ResidualPrior{make_variance(initial_value)};
}

auto make_marker_variance(
    MarkerVarianceLayout scope = MarkerVarianceLayout::per_marker,
    double initial_value = 1.0) -> MarkerVariance
{
    return MarkerVariance{scope, make_variance(initial_value)};
}

auto make_proportion_2() -> MixtureProportion
{
    return MixtureProportion{
        SimplexParameter{
            Eigen::VectorXd{{0.9, 0.1}},
            DirichletPrior{Eigen::VectorXd{{1.0, 1.0}}}},
        UpdatePolicy::fixed,
    };
}

auto make_proportion_4() -> MixtureProportion
{
    return MixtureProportion{
        SimplexParameter{
            Eigen::VectorXd{{0.91, 0.03, 0.03, 0.03}},
            DirichletPrior{Eigen::VectorXd{{1.0, 1.0, 1.0, 1.0}}}},
        UpdatePolicy::sampled,
    };
}

auto make_prior(std::vector<std::unique_ptr<GeneticPrior>> genetics)
    -> BayesPrior
{
    return BayesPrior(
        make_random_prior(0.25), std::move(genetics), make_residual_prior(1.5));
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

class FieldPathCollector final : public gelex::infra::FieldVisitor
{
   public:
    auto on(
        std::string_view name,
        Eigen::Ref<Eigen::VectorXf>,
        gelex::FieldFlag) -> void override
    {
        paths.push_back(make_path(name));
    }

    auto on(
        std::string_view name,
        Eigen::Ref<Eigen::VectorXd> value,
        gelex::FieldFlag) -> void override
    {
        const auto path = make_path(name);
        paths.push_back(path);
        if (path == mutate_path && value.size() > 0)
        {
            value(0) = mutated_value;
        }
    }

    auto on(
        std::string_view name,
        Eigen::Ref<Eigen::VectorXi>,
        gelex::FieldFlag) -> void override
    {
        paths.push_back(make_path(name));
    }

    auto on(std::string_view name, double& value, gelex::FieldFlag)
        -> void override
    {
        const auto path = make_path(name);
        paths.push_back(path);
        if (path == mutate_path)
        {
            value = mutated_value;
        }
    }

    auto on(std::string_view name, int&, gelex::FieldFlag) -> void override
    {
        paths.push_back(make_path(name));
    }

    auto count(std::string_view path) const -> std::size_t
    {
        return static_cast<std::size_t>(
            std::ranges::count(paths, std::string{path}));
    }

    std::vector<std::string> paths;
    std::string mutate_path;
    double mutated_value{9.0};

   private:
    auto enter(std::string_view name) -> void override
    {
        scopes_.emplace_back(name);
    }

    auto leave() -> void override { scopes_.pop_back(); }

    auto make_path(std::string_view name) const -> std::string
    {
        std::string path;
        for (const auto& scope : scopes_)
        {
            if (!path.empty())
            {
                path += '/';
            }
            path += scope;
        }
        if (!path.empty())
        {
            path += '/';
        }
        path += name;
        return path;
    }

    std::vector<std::string> scopes_;
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

TEST_CASE("BayesRecipe prior v2 constructs BayesStateV2", "[bayes_state]")
{
    constexpr std::array modes{GeneticMode::A};
    auto model = make_model(modes, true);

    BayesRecipeConfig config;
    config.modes = {GeneticMode::A};
    BayesRecipe recipe(BayesRecipePreset::RR, config);
    auto prior = recipe.make_prior_v2(model);

    BayesStateV2 state(model, prior);

    REQUIRE(state.fixed().coeffs.isApprox(Eigen::VectorXd::Zero(1)));
    REQUIRE(state.random().size() == 1);
    REQUIRE(state.random()[0].coeffs.size() == 2);
    REQUIRE(state.genetics().size() == 1);
    REQUIRE(state.genetic(GeneticMode::A) != nullptr);
    REQUIRE(state.genetic(GeneticMode::D) == nullptr);
    REQUIRE(state.residual().y_adj.isApprox(model.phenotype()));
    REQUIRE(state.residual().variance > 0.0);
}

TEST_CASE("BayesRecipe prior v2 constructs independent genetic blocks", "[bayes_state]")
{
    constexpr std::array modes{GeneticMode::A, GeneticMode::D};
    auto model = make_model(modes);

    BayesRecipeConfig config;
    config.modes = {GeneticMode::A, GeneticMode::D};
    BayesRecipe recipe(BayesRecipePreset::RR, config);
    auto prior = recipe.make_prior_v2(model);

    BayesStateV2 state(model, prior);

    REQUIRE(state.genetics().size() == 2);
    REQUIRE(state.genetic(GeneticMode::A) != nullptr);
    REQUIRE(state.genetic(GeneticMode::D) != nullptr);
    REQUIRE(
        state.genetic_block_for(GeneticMode::A)
        != state.genetic_block_for(GeneticMode::D));
    REQUIRE(
        std::holds_alternative<SingleGeneticBlockState>(
            *state.genetic_block_for(GeneticMode::A)));
    REQUIRE(
        std::holds_alternative<SingleGeneticBlockState>(
            *state.genetic_block_for(GeneticMode::D)));
}

TEST_CASE("BayesRecipe prior v2 constructs joint genetic blocks", "[bayes_state]")
{
    constexpr std::array modes{GeneticMode::A, GeneticMode::D};
    auto model = make_model(modes);

    BayesRecipeConfig config;
    config.modes = {GeneticMode::A, GeneticMode::D};
    BayesRecipe recipe(BayesRecipePreset::CD, config);
    auto prior = recipe.make_prior_v2(model);

    BayesStateV2 state(model, prior);

    REQUIRE(state.genetics().size() == 1);
    REQUIRE(state.genetic(GeneticMode::A) != nullptr);
    REQUIRE(state.genetic(GeneticMode::D) != nullptr);
    REQUIRE(
        state.genetic_block_for(GeneticMode::A)
        == state.genetic_block_for(GeneticMode::D));
    REQUIRE(
        std::holds_alternative<JointGeneticBlockState>(
            *state.genetic_block_for(GeneticMode::A)));
}

TEST_CASE(
    "BayesState initializes fixed random and residual state",
    "[bayes_state]")
{
    constexpr std::array modes{GeneticMode::A};
    auto model = make_model(modes, true);

    std::vector<std::unique_ptr<GeneticPrior>> genetics;
    genetics.push_back(
        std::make_unique<GaussianPrior>(
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

TEST_CASE(
    "FixedState RandomState ResidualState and GeneticState visit fields",
    "[bayes_state]")
{
    gelex::bayes::FixedState fixed{Eigen::VectorXd{{1.0, 2.0}}};
    FieldPathCollector fixed_visitor;
    fixed_visitor.mutate_path = "fixed/coeffs";

    fixed.visit(fixed_visitor);

    REQUIRE(fixed_visitor.count("fixed/coeffs") == 1);
    REQUIRE(fixed.coeffs(0) == fixed_visitor.mutated_value);

    gelex::bayes::RandomState random{Eigen::VectorXd{{3.0, 4.0}}, 0.5};
    FieldPathCollector random_visitor;
    random_visitor.mutate_path = "random/variance";

    random.visit(random_visitor);

    REQUIRE(random_visitor.count("random/coeffs") == 1);
    REQUIRE(random_visitor.count("random/variance") == 1);
    REQUIRE(random.variance == random_visitor.mutated_value);

    gelex::bayes::ResidualState residual{
        .y_adj = Eigen::VectorXd{{5.0, 6.0}},
        .variance = 0.75};
    FieldPathCollector residual_visitor;
    residual_visitor.mutate_path = "residual/y_adj";

    residual.visit(residual_visitor);

    REQUIRE(residual_visitor.count("residual/y_adj") == 1);
    REQUIRE(residual_visitor.count("residual/variance") == 1);
    REQUIRE(residual.y_adj(0) == residual_visitor.mutated_value);

    gelex::bayes::GeneticState genetic{GeneticMode::A, 2, 3};
    FieldPathCollector genetic_visitor;
    genetic_visitor.mutate_path = "genetic/coeffs";

    genetic.visit(genetic_visitor);

    REQUIRE(genetic_visitor.count("genetic/coeffs") == 1);
    REQUIRE(genetic_visitor.count("genetic/u") == 1);
    REQUIRE(genetic_visitor.count("genetic/variance") == 1);
    REQUIRE(genetic_visitor.count("genetic/heritability") == 1);
    REQUIRE(genetic.coeffs(0) == genetic_visitor.mutated_value);
}

TEST_CASE(
    "SingleGeneticBlockState and JointGeneticBlockState own genetic states",
    "[bayes_state]")
{
    STATIC_REQUIRE(std::variant_size_v<GeneticPriorBlockState> == 2);

    constexpr std::array modes{GeneticMode::A, GeneticMode::D};
    auto model = make_model(modes);

    SingleGaussianPrior single_prior{
        GeneticMode::A, make_marker_variance()};
    SingleGeneticBlockState single{
        *model.genetic(GeneticMode::A), single_prior};

    REQUIRE(single.mode() == GeneticMode::A);
    REQUIRE(single.contains(GeneticMode::A));
    REQUIRE(!single.contains(GeneticMode::D));
    REQUIRE(single.state().type == GeneticMode::A);
    REQUIRE(single.state().coeffs.size() == kNumMarkers);
    REQUIRE(single.state().u.size() == kNumIndividuals);
    REQUIRE(single.prior_state().query<SingleVarianceStateCap>() != nullptr);

    FieldPathCollector single_visitor;
    single_visitor.mutate_path = "single/genetic/coeffs";
    single.visit(single_visitor);

    REQUIRE(single_visitor.count("single/genetic/coeffs") == 1);
    REQUIRE(single_visitor.count("single/genetic/u") == 1);
    REQUIRE(single_visitor.count("single/prior_state/gaussian/variance") == 1);
    REQUIRE(single.state().coeffs(0) == single_visitor.mutated_value);

    JointGaussianMixturePrior joint_prior{
        std::array{
            make_marker_variance(MarkerVarianceLayout::per_marker, 0.5),
            make_marker_variance(MarkerVarianceLayout::shared, 0.75)},
        make_proportion_4()};
    JointGeneticBlockState joint{
        *model.genetic(GeneticMode::A),
        *model.genetic(GeneticMode::D),
        joint_prior};

    REQUIRE(joint.contains(GeneticMode::A));
    REQUIRE(joint.contains(GeneticMode::D));
    REQUIRE(joint.state(GeneticMode::A).type == GeneticMode::A);
    REQUIRE(joint.state(GeneticMode::D).type == GeneticMode::D);
    REQUIRE(joint.state(GeneticMode::A).coeffs.size() == kNumMarkers);
    REQUIRE(joint.state(GeneticMode::D).u.size() == kNumIndividuals);
    REQUIRE(joint.prior_state().query<JointVarianceStateCap>() != nullptr);

    FieldPathCollector joint_visitor;
    joint_visitor.mutate_path = "joint/D/genetic/variance";
    joint.visit(joint_visitor);

    REQUIRE(joint_visitor.count("joint/A/genetic/coeffs") == 1);
    REQUIRE(joint_visitor.count("joint/D/genetic/u") == 1);
    REQUIRE(
        joint_visitor.count(
            "joint/prior_state/joint_mixture_gaussian/A/variance")
        == 1);
    REQUIRE(joint.state(GeneticMode::D).variance == joint_visitor.mutated_value);
}

TEST_CASE(
    "BayesState creates genetic prior states by capability",
    "[bayes_state]")
{
    constexpr std::array modes{GeneticMode::A, GeneticMode::D};
    auto model = make_model(modes);

    SECTION("Gaussian")
    {
        std::vector<std::unique_ptr<GeneticPrior>> genetics;
        genetics.push_back(
            std::make_unique<GaussianPrior>(
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
        genetics.push_back(
            std::make_unique<SpikeSlabGaussianPrior>(
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
        REQUIRE(
            proportion.proportion()[0].count.isApprox(Eigen::VectorXi{{3, 0}}));
    }

    SECTION("ScaledMixture")
    {
        std::vector<std::unique_ptr<GeneticPrior>> genetics;
        genetics.push_back(
            std::make_unique<ScaledMixtureGaussianPrior>(
                GeneticMode::D,
                make_marker_variance(MarkerVarianceLayout::shared),
                Eigen::VectorXd{{0.0, 0.1, 1.0}},
                MixtureProportion{
                    SimplexParameter{
                        Eigen::VectorXd{{0.8, 0.1, 0.1}},
                        DirichletPrior{Eigen::VectorXd{{1.0, 1.0, 1.0}}}},
                    UpdatePolicy::fixed}));
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
    genetics.push_back(
        std::make_unique<JointMixtureGaussianPrior>(
            std::array{GeneticMode::A, GeneticMode::D},
            std::array{
                make_marker_variance(MarkerVarianceLayout::shared),
                make_marker_variance(MarkerVarianceLayout::shared)},
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
    REQUIRE(
        &state.genetics()[additive_block->genetic_indices()[0]] == additive);
    REQUIRE(
        &state.genetics()[additive_block->genetic_indices()[1]] == dominance);
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
    genetics.push_back(
        std::make_unique<GaussianPrior>(
            GeneticMode::D, make_marker_variance()));
    auto prior = make_prior(std::move(genetics));

    REQUIRE_THROWS_AS(BayesState(model, prior), GelexException);
}

TEST_CASE("BayesStateV2 creates single genetic blocks", "[bayes_state]")
{
    constexpr std::array modes{GeneticMode::A};
    auto model = make_model(modes, true);

    std::vector<GeneticPriorBlockV2> genetics;
    genetics.emplace_back(
        std::make_unique<SingleGaussianPrior>(
            GeneticMode::A, make_marker_variance()));
    BayesPriorV2 prior(
        make_random_prior(0.25), std::move(genetics), make_residual_prior(1.5));

    BayesStateV2 state(model, prior);

    REQUIRE(state.fixed().coeffs.isApprox(Eigen::VectorXd::Zero(1)));
    REQUIRE(state.random().size() == 1);
    REQUIRE(state.random()[0].variance == 0.25);
    REQUIRE(state.residual().variance == 1.5);
    REQUIRE(state.genetics().size() == 1);
    REQUIRE(state.genetic(GeneticMode::A) != nullptr);
    REQUIRE(state.genetic(GeneticMode::D) == nullptr);

    auto* block = state.genetic_block_for(GeneticMode::A);
    REQUIRE(block != nullptr);
    auto* single = std::get_if<SingleGeneticBlockState>(block);
    REQUIRE(single != nullptr);
    REQUIRE(single->state().coeffs.size() == kNumMarkers);
    REQUIRE(single->state().u.size() == kNumIndividuals);
    REQUIRE(single->prior_state().query<SingleVarianceStateCap>() != nullptr);

    state.genetic(GeneticMode::A)->variance = 3.0;
    state.residual().variance = 7.0;
    state.compute_heritability();

    REQUIRE(
        std::abs(
            state.genetic(GeneticMode::A)->heritability - 3.0 / 10.25)
        < 1e-12);
}

TEST_CASE("BayesStateV2 creates joint genetic blocks", "[bayes_state]")
{
    constexpr std::array modes{GeneticMode::A, GeneticMode::D};
    auto model = make_model(modes, true);

    std::vector<GeneticPriorBlockV2> genetics;
    genetics.emplace_back(
        std::make_unique<JointGaussianMixturePrior>(
            std::array{
                make_marker_variance(MarkerVarianceLayout::shared),
                make_marker_variance(MarkerVarianceLayout::shared)},
            make_proportion_4()));
    BayesPriorV2 prior(
        make_random_prior(2.0), std::move(genetics), make_residual_prior(10.0));

    BayesStateV2 state(model, prior);

    REQUIRE(state.genetics().size() == 1);
    REQUIRE(state.genetic(GeneticMode::A) != nullptr);
    REQUIRE(state.genetic(GeneticMode::D) != nullptr);
    REQUIRE(
        state.genetic_block_for(GeneticMode::A)
        == state.genetic_block_for(GeneticMode::D));

    auto* block = state.genetic_block_for(GeneticMode::A);
    REQUIRE(block != nullptr);
    auto* joint = std::get_if<JointGeneticBlockState>(block);
    REQUIRE(joint != nullptr);
    REQUIRE(joint->state(GeneticMode::A).coeffs.size() == kNumMarkers);
    REQUIRE(joint->state(GeneticMode::D).u.size() == kNumIndividuals);
    REQUIRE(joint->prior_state().query<JointVarianceStateCap>() != nullptr);

    state.genetic(GeneticMode::A)->variance = 3.0;
    state.genetic(GeneticMode::D)->variance = 5.0;
    state.compute_heritability();

    REQUIRE(
        std::abs(state.genetic(GeneticMode::A)->heritability - 0.15) < 1e-12);
    REQUIRE(
        std::abs(state.genetic(GeneticMode::D)->heritability - 0.25) < 1e-12);
}

TEST_CASE("BayesStateV2 visits fields", "[bayes_state]")
{
    constexpr std::array modes{GeneticMode::A};
    auto model = make_model(modes, true);

    std::vector<GeneticPriorBlockV2> genetics;
    genetics.emplace_back(
        std::make_unique<SingleGaussianPrior>(
            GeneticMode::A, make_marker_variance()));
    BayesPriorV2 prior(
        make_random_prior(0.25), std::move(genetics), make_residual_prior(1.5));

    BayesStateV2 state(model, prior);
    FieldPathCollector visitor;
    visitor.mutate_path = "state/random_0/variance";

    state.visit(visitor);

    REQUIRE(visitor.count("state/fixed/coeffs") == 1);
    REQUIRE(visitor.count("state/random_0/coeffs") == 1);
    REQUIRE(visitor.count("state/random_0/variance") == 1);
    REQUIRE(
        visitor.count("state/genetic_0/single/genetic/coeffs") == 1);
    REQUIRE(
        visitor.count(
            "state/genetic_0/single/prior_state/gaussian/variance")
        == 1);
    REQUIRE(visitor.count("state/residual/variance") == 1);
    REQUIRE(state.random()[0].variance == visitor.mutated_value);
}

TEST_CASE("BayesStateV2 rejects prior mode missing from model", "[bayes_state]")
{
    constexpr std::array modes{GeneticMode::A};
    auto model = make_model(modes);

    {
        std::vector<GeneticPriorBlockV2> genetics;
        genetics.emplace_back(
            std::make_unique<SingleGaussianPrior>(
                GeneticMode::D, make_marker_variance()));
        BayesPriorV2 prior(
            make_random_prior(), std::move(genetics), make_residual_prior());

        REQUIRE_THROWS_AS(BayesStateV2(model, prior), GelexException);
    }

    {
        std::vector<GeneticPriorBlockV2> genetics;
        genetics.emplace_back(
            std::make_unique<JointGaussianMixturePrior>(
                std::array{
                    make_marker_variance(MarkerVarianceLayout::shared),
                    make_marker_variance(MarkerVarianceLayout::shared)},
                make_proportion_4()));
        BayesPriorV2 prior(
            make_random_prior(), std::move(genetics), make_residual_prior());

        REQUIRE_THROWS_AS(BayesStateV2(model, prior), GelexException);
    }
}

TEST_CASE(
    "BayesState computes heritability from runtime variance",
    "[bayes_state]")
{
    constexpr std::array modes{GeneticMode::A, GeneticMode::D};
    auto model = make_model(modes, true);

    std::vector<std::unique_ptr<GeneticPrior>> genetics;
    genetics.push_back(
        std::make_unique<GaussianPrior>(
            GeneticMode::A, make_marker_variance()));
    genetics.push_back(
        std::make_unique<GaussianPrior>(
            GeneticMode::D, make_marker_variance()));
    auto prior = make_prior(std::move(genetics));

    BayesState state(model, prior);
    state.random()[0].variance = 2.0;
    state.genetic(GeneticMode::A)->variance = 3.0;
    state.genetic(GeneticMode::D)->variance = 5.0;
    state.residual().variance = 10.0;

    state.compute_heritability();

    REQUIRE(
        std::abs(state.genetic(GeneticMode::A)->heritability - 0.15) < 1e-12);
    REQUIRE(
        std::abs(state.genetic(GeneticMode::D)->heritability - 0.25) < 1e-12);
}

TEST_CASE(
    "BayesState visits sample records through state schema",
    "[bayes_state]")
{
    constexpr std::array modes{GeneticMode::A};
    auto model = make_model(modes, true);

    std::vector<std::unique_ptr<GeneticPrior>> genetics;
    genetics.push_back(
        std::make_unique<SpikeSlabGaussianPrior>(
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
    REQUIRE(
        sink.count("genetic_block/0/prior_state/proportion/0/assignment") == 1);
    REQUIRE(sink.count("residual/0/variance") == 1);
}

TEST_CASE(
    "BayesState visits joint prior state once for checkpoint",
    "[bayes_state]")
{
    constexpr std::array modes{GeneticMode::A, GeneticMode::D};
    auto model = make_model(modes);

    std::vector<std::unique_ptr<GeneticPrior>> genetics;
    genetics.push_back(
        std::make_unique<JointMixtureGaussianPrior>(
            std::array{GeneticMode::A, GeneticMode::D},
            std::array{
                make_marker_variance(MarkerVarianceLayout::shared),
                make_marker_variance(MarkerVarianceLayout::shared)},
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
    REQUIRE(
        sink.count("genetic_block/0/prior_state/proportion/0/assignment") == 1);
    REQUIRE(sink.count("genetic_block/0/prior_state/proportion/0/count") == 1);
    REQUIRE(sink.count("genetic_block/0/prior_state/proportion/0/value") == 1);
    REQUIRE(sink.count("genetic_block/1/prior_state/proportion/0/value") == 0);
}
