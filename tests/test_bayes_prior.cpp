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
#include <algorithm>
#include <memory>
#include <string>
#include <string_view>
#include <utility>
#include <variant>
#include <vector>

#include <Eigen/Core>
#include <catch2/catch_test_macros.hpp>

#include "gelex/exception.h"
#include "gelex/infra/field_flag.h"
#include "gelex/infra/field_visitor.h"
#include "gelex/model/bayes/capabilities.h"
#include "gelex/model/bayes/gaussian_prior.h"
#include "gelex/model/bayes/gaussian_prior_state.h"
#include "gelex/model/bayes/genetic_prior.h"
#include "gelex/model/bayes/genetic_priors/gaussian.h"
#include "gelex/model/bayes/prior.h"
#include "gelex/model/bayes/prior_parameters.h"
#include "gelex/model/bayes/prior_state.h"
#include "gelex/model/bayes/state_capabilities.h"
#include "gelex/types/genetic_effect_type.h"

using gelex::GelexException;
using gelex::GeneticMode;
using gelex::bayes::BayesPrior;
using gelex::bayes::BayesPriorV2;
using gelex::bayes::ComponentStateCap;
using gelex::bayes::DirichletPrior;
using gelex::bayes::GaussianPrior;
using gelex::bayes::GaussianState;
using gelex::bayes::GeneticPrior;
using gelex::bayes::GeneticPriorBlockV2;
using gelex::bayes::JointGeneticPrior;
using gelex::bayes::JointComponentStateCap;
using gelex::bayes::JointGaussianMixturePrior;
using gelex::bayes::JointGeneticPriorState;
using gelex::bayes::JointSharedMarkerVarianceCap;
using gelex::bayes::JointMixtureGaussianPrior;
using gelex::bayes::JointMixtureProportionCap;
using gelex::bayes::JointProportionStateCap;
using gelex::bayes::JointSharedVarianceStateCap;
using gelex::bayes::MarkerVariance;
using gelex::bayes::MarkerVarianceLayout;
using gelex::bayes::MixtureProportion;
using gelex::bayes::MixtureProportionCap;
using gelex::bayes::MultiplierCap;
using gelex::bayes::ProportionStateCap;
using gelex::bayes::RandomPrior;
using gelex::bayes::ResidualPrior;
using gelex::bayes::ScaledInvChiSqPrior;
using gelex::bayes::ScaledMixtureGaussianPrior;
using gelex::bayes::ScaledMixtureGaussianState;
using gelex::bayes::SharedMarkerVariance;
using gelex::bayes::SimplexParameter;
using gelex::bayes::SingleComponentStateCap;
using gelex::bayes::SingleGeneticPrior;
using gelex::bayes::SingleGeneticPriorState;
using gelex::bayes::SingleMixtureProportionCap;
using gelex::bayes::SingleMultiplierCap;
using gelex::bayes::SinglePerMarkerGaussianPrior;
using gelex::bayes::SinglePerMarkerSpikeSlabGaussianPrior;
using gelex::bayes::SinglePerMarkerVarianceCap;
using gelex::bayes::SinglePerMarkerVarianceStateCap;
using gelex::bayes::SingleProportionStateCap;
using gelex::bayes::SingleScaledMixtureGaussianPrior;
using gelex::bayes::SingleSharedGaussianPrior;
using gelex::bayes::SingleSharedMarkerVarianceCap;
using gelex::bayes::SingleSharedSpikeSlabGaussianPrior;
using gelex::bayes::SingleSharedVarianceStateCap;
using gelex::bayes::SpikeSlabGaussianPrior;
using gelex::bayes::SpikeSlabGaussianState;
using gelex::bayes::UpdatePolicy;
using gelex::bayes::VarianceParameter;
using gelex::bayes::JointSharedMarkerVariance;
using gelex::bayes::JointMixtureGaussianState;
using gelex::bayes::PerMarkerVariance;

namespace
{

auto make_variance(double init = 1.0) -> VarianceParameter
{
    return VarianceParameter(init, ScaledInvChiSqPrior{4.0, 1.0});
}

auto make_random_prior(double init = 1.0) -> RandomPrior
{
    return RandomPrior{make_variance(init)};
}

auto make_residual_prior(double init = 1.0) -> ResidualPrior
{
    return ResidualPrior{make_variance(init)};
}

auto make_marker_variance() -> MarkerVariance
{
    return MarkerVariance{
        MarkerVarianceLayout::per_marker,
        make_variance(),
    };
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

auto make_multiplier_2() -> Eigen::VectorXd
{
    return Eigen::VectorXd{{0.0, 1.0}};
}

auto make_gaussian(GeneticMode mode = GeneticMode::A)
    -> std::unique_ptr<GeneticPrior>
{
    return std::make_unique<GaussianPrior>(mode, make_marker_variance());
}

auto make_spike_slab(GeneticMode mode = GeneticMode::A)
    -> std::unique_ptr<GeneticPrior>
{
    return std::make_unique<SpikeSlabGaussianPrior>(
        mode, make_marker_variance(), make_proportion_2());
}

auto make_joint() -> std::unique_ptr<GeneticPrior>
{
    return std::make_unique<JointMixtureGaussianPrior>(
        std::array<GeneticMode, 2>{GeneticMode::A, GeneticMode::D},
        std::array<MarkerVariance, 2>{
            make_marker_variance(), make_marker_variance()},
        make_proportion_2());
}

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

    auto leave() -> void override
    {
        scopes_.pop_back();
    }

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

TEST_CASE("BayesPrior accessors expose construction arguments", "[bayes_prior]")
{
    std::vector<std::unique_ptr<GeneticPrior>> genetics;
    genetics.push_back(make_gaussian(GeneticMode::A));
    genetics.push_back(make_spike_slab(GeneticMode::D));

    BayesPrior prior(
        make_random_prior(2.0), std::move(genetics), make_residual_prior(5.0));

    REQUIRE(prior.random().initial_value() == 2.0);
    REQUIRE(prior.residual().initial_value() == 5.0);

    const auto range = prior.genetics();
    REQUIRE(range.size() == 2);
    REQUIRE_FALSE(range.empty());

    auto it = range.begin();
    REQUIRE((*it).contains(GeneticMode::A));
    ++it;
    REQUIRE((*it).contains(GeneticMode::D));
    ++it;
    REQUIRE(it == range.end());
}

TEST_CASE("BayesPriorV2 accessors expose genetic blocks", "[bayes_prior]")
{
    std::vector<GeneticPriorBlockV2> genetics;
    genetics.emplace_back(
        std::make_unique<SingleSharedGaussianPrior>(
            GeneticMode::A, SharedMarkerVariance{make_variance()}));
    genetics.emplace_back(
        std::make_unique<SinglePerMarkerSpikeSlabGaussianPrior>(
            GeneticMode::D,
            PerMarkerVariance{make_variance()},
            make_proportion_2()));

    BayesPriorV2 prior(
        make_random_prior(2.0), std::move(genetics), make_residual_prior(5.0));

    REQUIRE(prior.random().initial_value() == 2.0);
    REQUIRE(prior.residual().initial_value() == 5.0);

    const auto blocks = prior.genetics();
    REQUIRE(blocks.size() == 2);
    const auto* additive
        = std::get_if<std::unique_ptr<SingleGeneticPrior>>(&blocks[0]);
    const auto* dominance
        = std::get_if<std::unique_ptr<SingleGeneticPrior>>(&blocks[1]);
    REQUIRE(additive != nullptr);
    REQUIRE((*additive)->mode() == GeneticMode::A);
    REQUIRE(dominance != nullptr);
    REQUIRE((*dominance)->mode() == GeneticMode::D);
}

TEST_CASE("BayesPriorV2 accepts a joint genetic block", "[bayes_prior]")
{
    std::vector<GeneticPriorBlockV2> genetics;
    genetics.emplace_back(
        std::make_unique<JointGaussianMixturePrior>(
            JointSharedMarkerVariance{std::array{
                SharedMarkerVariance{make_variance()},
                SharedMarkerVariance{make_variance()}}},
            make_proportion_2()));

    BayesPriorV2 prior(
        make_random_prior(), std::move(genetics), make_residual_prior());

    const auto blocks = prior.genetics();
    REQUIRE(blocks.size() == 1);
    REQUIRE(
        std::holds_alternative<std::unique_ptr<JointGeneticPrior>>(blocks[0]));
}

TEST_CASE("BayesPriorV2 visits prior fields", "[bayes_prior]")
{
    std::vector<GeneticPriorBlockV2> genetics;
    genetics.emplace_back(
        std::make_unique<SingleSharedGaussianPrior>(
            GeneticMode::A, SharedMarkerVariance{make_variance()}));
    genetics.emplace_back(
        std::make_unique<SinglePerMarkerSpikeSlabGaussianPrior>(
            GeneticMode::D,
            PerMarkerVariance{make_variance()},
            make_proportion_2()));
    BayesPriorV2 prior(
        make_random_prior(), std::move(genetics), make_residual_prior());
    FieldPathCollector visitor;

    prior.visit(visitor);

    REQUIRE(visitor.count("prior/random/variance/initial_value") == 1);
    REQUIRE(visitor.count("prior/genetic/0/shared_gaussian/mode") == 1);
    REQUIRE(
        visitor.count(
            "prior/genetic/0/shared_gaussian/shared_marker_variance/variance/"
            "initial_value")
        == 1);
    REQUIRE(
        visitor.count("prior/genetic/1/per_marker_spike_slab_gaussian/mode")
        == 1);
    REQUIRE(visitor.count("prior/residual/variance/initial_value") == 1);
}

TEST_CASE("BayesPriorV2 rejects null genetic blocks", "[bayes_prior]")
{
    std::vector<GeneticPriorBlockV2> genetics;

    SECTION("single")
    {
        genetics.emplace_back(std::unique_ptr<SingleGeneticPrior>{});
        REQUIRE_THROWS_AS(
            BayesPriorV2(
                make_random_prior(),
                std::move(genetics),
                make_residual_prior()),
            GelexException);
    }

    SECTION("joint")
    {
        genetics.emplace_back(std::unique_ptr<JointGeneticPrior>{});
        REQUIRE_THROWS_AS(
            BayesPriorV2(
                make_random_prior(),
                std::move(genetics),
                make_residual_prior()),
            GelexException);
    }
}

TEST_CASE("BayesPriorV2 rejects duplicate genetic modes", "[bayes_prior]")
{
    SECTION("single single conflict")
    {
        std::vector<GeneticPriorBlockV2> genetics;
        genetics.emplace_back(
            std::make_unique<SingleSharedGaussianPrior>(
                GeneticMode::A, SharedMarkerVariance{make_variance()}));
        genetics.emplace_back(
            std::make_unique<SinglePerMarkerSpikeSlabGaussianPrior>(
                GeneticMode::A,
                PerMarkerVariance{make_variance()},
                make_proportion_2()));

        REQUIRE_THROWS_AS(
            BayesPriorV2(
                make_random_prior(),
                std::move(genetics),
                make_residual_prior()),
            GelexException);
    }

    SECTION("single joint conflict")
    {
        std::vector<GeneticPriorBlockV2> genetics;
        genetics.emplace_back(
            std::make_unique<SingleSharedGaussianPrior>(
                GeneticMode::D, SharedMarkerVariance{make_variance()}));
        genetics.emplace_back(
            std::make_unique<JointGaussianMixturePrior>(
                JointSharedMarkerVariance{std::array{
                    SharedMarkerVariance{make_variance()},
                    SharedMarkerVariance{make_variance()}}},
                make_proportion_2()));

        REQUIRE_THROWS_AS(
            BayesPriorV2(
                make_random_prior(),
                std::move(genetics),
                make_residual_prior()),
            GelexException);
    }
}

TEST_CASE("GeneticPrior capabilities compose prior data axes", "[bayes_prior]")
{
    GaussianPrior gaussian(GeneticMode::A, make_marker_variance());
    REQUIRE(gaussian.variance().size() == 1);
    REQUIRE(gaussian.query<MixtureProportionCap>() == nullptr);
    REQUIRE(gaussian.query<MultiplierCap>() == nullptr);

    SpikeSlabGaussianPrior spike_slab(
        GeneticMode::A, make_marker_variance(), make_proportion_2());
    const auto* spike_slab_proportion
        = spike_slab.query<MixtureProportionCap>();
    REQUIRE(spike_slab.variance().size() == 1);
    REQUIRE(spike_slab_proportion != nullptr);
    REQUIRE(spike_slab_proportion->proportion().size() == 1);
    REQUIRE(spike_slab.query<MultiplierCap>() == nullptr);

    ScaledMixtureGaussianPrior scaled_mixture(
        GeneticMode::A,
        make_marker_variance(),
        make_multiplier_2(),
        make_proportion_2());
    const auto* scaled_mixture_multiplier
        = scaled_mixture.query<MultiplierCap>();
    REQUIRE(scaled_mixture.variance().size() == 1);
    REQUIRE(scaled_mixture.query<MixtureProportionCap>() != nullptr);
    REQUIRE(scaled_mixture_multiplier != nullptr);
    REQUIRE(scaled_mixture_multiplier->multiplier().size() == 1);

    JointMixtureGaussianPrior joint(
        std::array<GeneticMode, 2>{GeneticMode::A, GeneticMode::D},
        std::array<MarkerVariance, 2>{
            make_marker_variance(), make_marker_variance()},
        make_proportion_2());
    const auto* joint_proportion = joint.query<MixtureProportionCap>();
    REQUIRE(joint.variance().size() == 2);
    REQUIRE(joint_proportion != nullptr);
    REQUIRE(joint_proportion->proportion().size() == 1);
    REQUIRE(joint.query<MultiplierCap>() == nullptr);
}

TEST_CASE(
    "Single and joint genetic prior capabilities expose distinct shapes",
    "[bayes_prior]")
{
    SingleSharedGaussianPrior shared_gaussian(
        GeneticMode::A, SharedMarkerVariance{make_variance()});
    const auto* shared_variance
        = shared_gaussian.query<SingleSharedMarkerVarianceCap>();
    REQUIRE(shared_variance != nullptr);
    REQUIRE(shared_variance->variance().parameter().initial_value() == 1.0);
    REQUIRE(shared_gaussian.query<SinglePerMarkerVarianceCap>() == nullptr);
    REQUIRE(shared_gaussian.query<SingleMixtureProportionCap>() == nullptr);
    REQUIRE(shared_gaussian.query<SingleMultiplierCap>() == nullptr);
    REQUIRE(shared_gaussian.query<JointSharedMarkerVarianceCap>() == nullptr);

    SinglePerMarkerGaussianPrior per_marker_gaussian(
        GeneticMode::A, PerMarkerVariance{make_variance()});
    const auto* per_marker_variance
        = per_marker_gaussian.query<SinglePerMarkerVarianceCap>();
    REQUIRE(per_marker_variance != nullptr);
    REQUIRE(per_marker_variance->variance().parameter().initial_value() == 1.0);
    REQUIRE(
        per_marker_gaussian.query<SingleSharedMarkerVarianceCap>() == nullptr);

    SingleScaledMixtureGaussianPrior scaled_mixture(
        GeneticMode::D,
        SharedMarkerVariance{make_variance()},
        make_multiplier_2(),
        make_proportion_2());
    const auto* single_variance
        = scaled_mixture.query<SingleSharedMarkerVarianceCap>();
    const auto* single_proportion
        = scaled_mixture.query<SingleMixtureProportionCap>();
    const auto* single_multiplier
        = scaled_mixture.query<SingleMultiplierCap>();
    REQUIRE(single_variance != nullptr);
    REQUIRE(single_proportion != nullptr);
    REQUIRE(single_proportion->proportion().size() == 2);
    REQUIRE(single_multiplier != nullptr);
    REQUIRE(single_multiplier->multiplier().isApprox(make_multiplier_2()));
    REQUIRE(scaled_mixture.query<JointMixtureProportionCap>() == nullptr);

    JointGaussianMixturePrior joint(
        JointSharedMarkerVariance{std::array{
            SharedMarkerVariance{make_variance()},
            SharedMarkerVariance{make_variance()}}},
        make_proportion_2());
    const auto* joint_variance = joint.query<JointSharedMarkerVarianceCap>();
    const auto* joint_proportion = joint.query<JointMixtureProportionCap>();
    REQUIRE(joint_variance != nullptr);
    REQUIRE(
        &joint_variance->variance(GeneticMode::A)
        != &joint_variance->variance(GeneticMode::D));
    REQUIRE(joint_proportion != nullptr);
    REQUIRE(joint_proportion->proportion().size() == 2);
    REQUIRE(joint.query<SingleSharedMarkerVarianceCap>() == nullptr);
    REQUIRE(joint.query<SinglePerMarkerVarianceCap>() == nullptr);
}

TEST_CASE(
    "Single and joint genetic prior states expose distinct shapes",
    "[bayes_prior]")
{
    constexpr Eigen::Index kNumMarkers = 3;
    constexpr Eigen::Index kNumIndividuals = 5;

    SingleSharedGaussianPrior shared_gaussian(
        GeneticMode::A, SharedMarkerVariance{make_variance()});
    auto shared_state
        = shared_gaussian.make_state(kNumMarkers, kNumIndividuals);
    auto& shared_variance
        = shared_state->require<SingleSharedVarianceStateCap>();
    REQUIRE(shared_state->query<SingleGeneticPriorState>() != nullptr);
    REQUIRE(shared_state->query<JointSharedVarianceStateCap>() == nullptr);
    REQUIRE(shared_variance.variance() == 1.0);
    REQUIRE(shared_state->query<SinglePerMarkerVarianceStateCap>() == nullptr);
    REQUIRE(shared_state->query<SingleProportionStateCap>() == nullptr);

    SinglePerMarkerGaussianPrior per_marker_gaussian(
        GeneticMode::A, PerMarkerVariance{make_variance()});
    auto per_marker_state
        = per_marker_gaussian.make_state(kNumMarkers, kNumIndividuals);
    auto& per_marker_variance
        = per_marker_state->require<SinglePerMarkerVarianceStateCap>();
    REQUIRE(per_marker_variance.variance().isApprox(
        Eigen::VectorXd::Constant(kNumMarkers, 1.0)));
    REQUIRE(per_marker_state->query<SingleSharedVarianceStateCap>() == nullptr);

    SingleScaledMixtureGaussianPrior scaled_mixture(
        GeneticMode::D,
        SharedMarkerVariance{make_variance()},
        make_multiplier_2(),
        make_proportion_2());
    auto scaled_mixture_state
        = scaled_mixture.make_state(kNumMarkers, kNumIndividuals);
    auto& single_component
        = scaled_mixture_state->require<SingleComponentStateCap>();
    auto& single_proportion
        = scaled_mixture_state->require<SingleProportionStateCap>();
    REQUIRE(scaled_mixture_state->query<JointGeneticPriorState>() == nullptr);
    REQUIRE(single_component.component().gebv.size() == 1);
    REQUIRE(single_component.component().gebv[0].size() == kNumIndividuals);
    REQUIRE(single_component.component().gebv_var.size() == 1);
    REQUIRE(single_proportion.proportion().assignment.size() == kNumMarkers);
    REQUIRE(single_proportion.proportion().count.isApprox(
        Eigen::VectorXi{{3, 0}}));

    JointGaussianMixturePrior joint(
        JointSharedMarkerVariance{std::array{
            SharedMarkerVariance{make_variance()},
            SharedMarkerVariance{make_variance()}}},
        make_proportion_2());
    auto joint_state = joint.make_state(4, kNumIndividuals);
    auto& joint_variance = joint_state->require<JointSharedVarianceStateCap>();
    auto& joint_component = joint_state->require<JointComponentStateCap>();
    auto& joint_proportion = joint_state->require<JointProportionStateCap>();
    REQUIRE(joint_state->query<JointGeneticPriorState>() != nullptr);
    REQUIRE(joint_state->query<SingleSharedVarianceStateCap>() == nullptr);
    REQUIRE(joint_state->query<SinglePerMarkerVarianceStateCap>() == nullptr);
    REQUIRE(joint_state->query<ProportionStateCap>() == nullptr);
    REQUIRE(
        &joint_variance.variance(GeneticMode::A)
        != &joint_variance.variance(GeneticMode::D));
    REQUIRE(joint_variance.variance(GeneticMode::A) == 1.0);
    REQUIRE(joint_variance.variance(GeneticMode::D) == 1.0);
    REQUIRE(joint_component.component().gebv.size() == 1);
    REQUIRE(joint_component.component().gebv[0].size() == kNumIndividuals);
    REQUIRE(joint_proportion.proportion().assignment.size() == 4);
}

TEST_CASE("Single and joint genetic prior states visit fields", "[bayes_prior]")
{
    constexpr Eigen::Index kNumMarkers = 3;
    constexpr Eigen::Index kNumIndividuals = 5;

    SingleSharedGaussianPrior gaussian(
        GeneticMode::A, SharedMarkerVariance{make_variance()});
    auto gaussian_state = gaussian.make_state(kNumMarkers, kNumIndividuals);
    FieldPathCollector gaussian_visitor;
    gaussian_visitor.mutate_path = "shared_gaussian/variance";

    gaussian_state->visit(gaussian_visitor);

    REQUIRE(gaussian_visitor.count("shared_gaussian/variance") == 1);
    REQUIRE(
        gaussian_state->require<SingleSharedVarianceStateCap>().variance()
        == gaussian_visitor.mutated_value);

    SingleScaledMixtureGaussianPrior scaled_mixture(
        GeneticMode::D,
        SharedMarkerVariance{make_variance()},
        make_multiplier_2(),
        make_proportion_2());
    auto scaled_mixture_state
        = scaled_mixture.make_state(kNumMarkers, kNumIndividuals);
    FieldPathCollector scaled_mixture_visitor;

    scaled_mixture_state->visit(scaled_mixture_visitor);

    REQUIRE(
        scaled_mixture_visitor.count("scaled_mixture_gaussian/variance")
        == 1);
    REQUIRE(
        scaled_mixture_visitor.count(
            "scaled_mixture_gaussian/component/gebv_var")
        == 1);
    REQUIRE(
        scaled_mixture_visitor.count(
            "scaled_mixture_gaussian/component/gebv_0")
        == 1);
    REQUIRE(
        scaled_mixture_visitor.count(
            "scaled_mixture_gaussian/proportion/assignment")
        == 1);
    REQUIRE(
        scaled_mixture_visitor.count(
            "scaled_mixture_gaussian/proportion/update")
        == 1);

    JointGaussianMixturePrior joint(
        JointSharedMarkerVariance{std::array{
            SharedMarkerVariance{make_variance()},
            SharedMarkerVariance{make_variance()}}},
        make_proportion_2());
    auto joint_state = joint.make_state(kNumMarkers, kNumIndividuals);
    FieldPathCollector joint_visitor;

    joint_state->visit(joint_visitor);

    REQUIRE(joint_visitor.count("joint_mixture_gaussian/A/variance") == 1);
    REQUIRE(joint_visitor.count("joint_mixture_gaussian/D/variance") == 1);
    REQUIRE(
        joint_visitor.count("joint_mixture_gaussian/component/gebv_var") == 1);
    REQUIRE(
        joint_visitor.count("joint_mixture_gaussian/proportion/value") == 1);
}

TEST_CASE(
    "GeneticPrior::make_state builds capability-composed state",
    "[bayes_prior]")
{
    GaussianPrior gaussian(GeneticMode::A, make_marker_variance());
    constexpr Eigen::Index kSingleNumMarkers = 3;
    constexpr Eigen::Index kNumIndividuals = 5;

    auto gaussian_state
        = gaussian.make_state(kSingleNumMarkers, kNumIndividuals);
    auto& gaussian_variance = gaussian_state->require<GaussianState>();
    REQUIRE(gaussian_variance.variance().size() == 1);
    REQUIRE(gaussian_variance.variance()[0].size() == 3);
    REQUIRE(gaussian_variance.variance()[0].isApprox(
        Eigen::VectorXd::Constant(3, 1.0)));
    REQUIRE(gaussian_state->query<ProportionStateCap>() == nullptr);

    SpikeSlabGaussianPrior spike_slab(
        GeneticMode::A, make_marker_variance(), make_proportion_2());
    auto spike_slab_state
        = spike_slab.make_state(kSingleNumMarkers, kNumIndividuals);
    auto& spike_slab_variance
        = spike_slab_state->require<SpikeSlabGaussianState>();
    auto& spike_slab_proportion
        = spike_slab_state->require<ProportionStateCap>();
    REQUIRE(spike_slab_variance.variance().size() == 1);
    REQUIRE(spike_slab_proportion.proportion().size() == 1);
    REQUIRE(spike_slab_proportion.proportion()[0].assignment.size() == 3);
    REQUIRE(spike_slab_proportion.proportion()[0].assignment.isZero());
    REQUIRE(spike_slab_proportion.proportion()[0].count.isApprox(
        Eigen::VectorXi{{3, 0}}));
    REQUIRE(spike_slab_proportion.proportion()[0].value.isApprox(
        Eigen::VectorXd{{0.9, 0.1}}));

    ScaledMixtureGaussianPrior scaled_mixture(
        GeneticMode::A,
        make_marker_variance(),
        make_multiplier_2(),
        make_proportion_2());
    auto scaled_mixture_state
        = scaled_mixture.make_state(kSingleNumMarkers, kNumIndividuals);
    auto& scaled_mixture_component
        = scaled_mixture_state->require<ComponentStateCap>();
    auto& scaled_mixture_proportion
        = scaled_mixture_state->require<ProportionStateCap>();
    REQUIRE(scaled_mixture_component.component().size() == 1);
    REQUIRE(scaled_mixture_component.component()[0].gebv.size() == 1);
    REQUIRE(
        scaled_mixture_component.component()[0].gebv[0].size()
        == kNumIndividuals);
    REQUIRE(scaled_mixture_component.component()[0].gebv_var.size() == 1);
    REQUIRE(scaled_mixture_proportion.proportion().size() == 1);
    REQUIRE(scaled_mixture_proportion.proportion()[0].assignment.size() == 3);
    auto& scaled_mixture_variance
        = scaled_mixture_state->require<ScaledMixtureGaussianState>();
    REQUIRE(scaled_mixture_variance.variance().size() == 1);
    REQUIRE(scaled_mixture_variance.variance()[0].size() == kSingleNumMarkers);

    JointMixtureGaussianPrior joint(
        std::array<GeneticMode, 2>{GeneticMode::A, GeneticMode::D},
        std::array<MarkerVariance, 2>{
            make_marker_variance(), make_marker_variance()},
        make_proportion_2());

    auto joint_state = joint.make_state(4, kNumIndividuals);
    auto& joint_variance = joint_state->require<JointMixtureGaussianState>();
    auto& joint_component = joint_state->require<ComponentStateCap>();
    auto& joint_proportion = joint_state->require<ProportionStateCap>();
    REQUIRE(joint_variance.variance().size() == 2);
    REQUIRE(joint_variance.variance()[0].size() == 4);
    REQUIRE(joint_variance.variance()[1].size() == 4);
    REQUIRE(joint_component.component().size() == 1);
    REQUIRE(joint_component.component()[0].gebv.size() == 1);
    REQUIRE(joint_component.component()[0].gebv[0].size() == kNumIndividuals);
    REQUIRE(joint_component.component()[0].gebv_var.size() == 1);
    REQUIRE(joint_proportion.proportion().size() == 1);
    REQUIRE(joint_proportion.proportion()[0].assignment.size() == 4);
}

TEST_CASE("BayesPrior::genetics range supports for-each", "[bayes_prior]")
{
    std::vector<std::unique_ptr<GeneticPrior>> genetics;
    genetics.push_back(make_gaussian(GeneticMode::A));
    genetics.push_back(make_spike_slab(GeneticMode::D));

    BayesPrior prior(
        make_random_prior(), std::move(genetics), make_residual_prior());

    std::vector<GeneticMode> visited;
    for (const auto& block : prior.genetics())
    {
        for (const auto m : block.modes())
        {
            visited.push_back(m);
        }
    }
    REQUIRE(
        visited == std::vector<GeneticMode>{GeneticMode::A, GeneticMode::D});
}

TEST_CASE(
    "BayesPrior::genetics empty for prior with no blocks",
    "[bayes_prior]")
{
    BayesPrior prior(make_random_prior(), {}, make_residual_prior());

    const auto range = prior.genetics();
    REQUIRE(range.empty());
    REQUIRE(range.size() == 0);
    REQUIRE(range.begin() == range.end());
}

TEST_CASE("BayesPrior rejects null genetic block", "[bayes_prior]")
{
    std::vector<std::unique_ptr<GeneticPrior>> genetics;
    genetics.push_back(nullptr);

    REQUIRE_THROWS_AS(
        BayesPrior(
            make_random_prior(), std::move(genetics), make_residual_prior()),
        GelexException);
}

TEST_CASE("BayesPrior rejects duplicate mode across blocks", "[bayes_prior]")
{
    std::vector<std::unique_ptr<GeneticPrior>> genetics;
    genetics.push_back(make_gaussian(GeneticMode::A));
    genetics.push_back(make_gaussian(GeneticMode::A));

    REQUIRE_THROWS_AS(
        BayesPrior(
            make_random_prior(), std::move(genetics), make_residual_prior()),
        GelexException);
}

TEST_CASE(
    "VarianceParameter rejects invalid initial_value",
    "[variance_parameter]")
{
    SECTION("initial_value <= 0")
    {
        REQUIRE_THROWS_AS(
            VarianceParameter(0.0, ScaledInvChiSqPrior{4.0, 1.0}),
            GelexException);
        REQUIRE_THROWS_AS(
            VarianceParameter(-1.0, ScaledInvChiSqPrior{4.0, 1.0}),
            GelexException);
    }
}

TEST_CASE(
    "VarianceParameter accepts positive initial_value",
    "[variance_parameter]")
{
    const VarianceParameter parameter{2.0, ScaledInvChiSqPrior{4.0, 1.0}};

    REQUIRE(parameter.initial_value() == 2.0);
}

TEST_CASE(
    "SimplexParameter rejects invalid initial_value",
    "[simplex_parameter]")
{
    SECTION("initial_value size < 2")
    {
        REQUIRE_THROWS_AS(
            SimplexParameter(
                Eigen::VectorXd{{1.0}},
                DirichletPrior{Eigen::VectorXd{{1.0, 1.0}}}),
            GelexException);
    }

    SECTION("initial_value entries must be positive")
    {
        REQUIRE_THROWS_AS(
            SimplexParameter(
                Eigen::VectorXd{{1.0, 0.0}},
                DirichletPrior{Eigen::VectorXd{{1.0, 1.0}}}),
            GelexException);
    }

    SECTION("initial_value entries must sum to one")
    {
        REQUIRE_THROWS_AS(
            SimplexParameter(
                Eigen::VectorXd{{0.7, 0.2}},
                DirichletPrior{Eigen::VectorXd{{1.0, 1.0}}}),
            GelexException);
    }

    SECTION("prior size must match initial_value")
    {
        REQUIRE_THROWS_AS(
            SimplexParameter(
                Eigen::VectorXd{{0.5, 0.5}},
                DirichletPrior{Eigen::VectorXd{{1.0, 1.0, 1.0}}}),
            GelexException);
    }
}

TEST_CASE("DirichletPrior rejects invalid concentration", "[dirichlet_prior]")
{
    SECTION("concentration entries must be positive")
    {
        REQUIRE_THROWS_AS(
            DirichletPrior(Eigen::VectorXd{{1.0, 0.0}}), GelexException);
    }

    SECTION("concentration size < 2")
    {
        REQUIRE_THROWS_AS(
            DirichletPrior(Eigen::VectorXd{{1.0}}), GelexException);
    }
}

TEST_CASE("ScaledInvChiSqPrior rejects invalid inputs", "[chisq_prior]")
{
    SECTION("degrees_of_freedom <= 0")
    {
        REQUIRE_THROWS_AS(ScaledInvChiSqPrior(0.0, 1.0), GelexException);
    }
    SECTION("scale <= 0 for proper prior")
    {
        REQUIRE_THROWS_AS(ScaledInvChiSqPrior(4.0, 0.0), GelexException);
    }
}

TEST_CASE(
    "ScaledInvChiSqPrior accepts flat prior sentinel {-2, 0}",
    "[chisq_prior]")
{
    SECTION("flat sentinel is accepted")
    {
        REQUIRE_NOTHROW(ScaledInvChiSqPrior(-2.0, 0.0));
    }
}

TEST_CASE(
    "ScaledInvChiSqPrior rejects near-flat but invalid prior combos",
    "[chisq_prior]")
{
    SECTION("{-2, 1}: df matches sentinel but scale != 0")
    {
        REQUIRE_THROWS_AS(ScaledInvChiSqPrior(-2.0, 1.0), GelexException);
    }
    SECTION("{4, 0}: scale matches sentinel but df != -2")
    {
        REQUIRE_THROWS_AS(ScaledInvChiSqPrior(4.0, 0.0), GelexException);
    }
    SECTION("{-1, 0}: negative df but not sentinel, scale == 0")
    {
        REQUIRE_THROWS_AS(ScaledInvChiSqPrior(-1.0, 0.0), GelexException);
    }
}
