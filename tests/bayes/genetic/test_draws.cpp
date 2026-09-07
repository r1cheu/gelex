// Copyright 2026 RuLei Chen
// SPDX-License-Identifier: Apache-2.0

#include <Eigen/Core>
#include <array>
#include <catch2/catch_test_macros.hpp>
#include <cstdint>
#include <string>

#include "gelex/bayes/genetic/construction.h"
#include "gelex/bayes/genetic/gaussian.h"
#include "gelex/bayes/genetic/joint_spike_slab.h"
#include "gelex/bayes/genetic/parameter.h"
#include "gelex/bayes/genetic/scaled_mixture.h"
#include "gelex/bayes/genetic/spike_slab.h"
#include "gelex/bayes/genetic/types.h"
#include "gelex/bayes/mode_values.h"
#include "gelex/bayes/prior.h"
#include "gelex/bayes/recipe.h"
#include "gelex/bayes/spec.h"
#include "gelex/bayes/stats/dirichlet_log_kernel.h"
#include "gelex/bayes/stats/scaled_inv_chi2_log_kernel.h"
#include "gelex/bayes/variance/budget.h"
#include "gelex/genetic_mode.h"
#include "gelex/io/binary_reader.h"
#include "gelex/io/binary_writer.h"

#include "bayes/bayes_model_fixture.h"
#include "file_fixture.h"

namespace
{

template <typename Spec>
auto check_mode_draws(const Spec& spec) -> void
{
    constexpr auto modes = gelex::GeneticModeSet{gelex::GeneticMode::A};
    const auto model = gelex::test::make_random_effect_model(modes);
    const auto full_prior = gelex::make_prior(
        gelex::BayesRecipe<modes, gelex::HomogeneousModeValues<modes, Spec>>{
            gelex::HomogeneousModeValues<modes, Spec>{spec},
            gelex::VarianceBudget{{.additive = 0.4, .random = 0.1}}},
        model);
    const auto& prior
        = full_prior.genetic().template get<gelex::GeneticMode::A>();
    gelex::test::FileFixture fixture;
    const auto path = (fixture.get_test_dir() / "mode.draws").string();
    const gelex::GeneticDimensions dimensions{.individual = 3, .marker = 2};
    auto state = gelex::detail::make_state(prior, dimensions);
    {
        gelex::BinaryWriter writer{path};
        auto draws = gelex::detail::make_draws(
            prior, writer, "genetic/A", 2, dimensions);
        if constexpr (requires { state.assignments(); })
        {
            static_cast<void>(state.transition(0, 1.25, 1));
        }
        else
        {
            state.transition(0, 1.25);
        }
        draws.append(state);
        if constexpr (requires { state.assignments(); })
        {
            static_cast<void>(state.transition(0, 0.0, 0));
            static_cast<void>(state.transition(1, -2.5, 1));
        }
        else
        {
            state.transition(0, 0.0);
            state.transition(1, -2.5);
        }
        draws.append(state);
        writer.close();
    }
    const gelex::BinaryReader reader{path};
    REQUIRE(reader.to_map<float>("genetic/A/coefficients")
                .isApprox(Eigen::MatrixXf{{1.25F, 0.0F}, {0.0F, -2.5F}}));
    if constexpr (requires { state.assignments(); })
    {
        REQUIRE(reader.to_map<std::uint8_t>("genetic/A/assignment")
                    .cast<double>()
                    .isApprox(Eigen::MatrixXd{{1, 0}, {0, 1}}));
    }
    if constexpr (requires { state.variance().size(); })
    {
        REQUIRE(reader.to_map<float>("genetic/A/variance")
                    .isApprox(
                        Eigen::MatrixXf::Constant(
                            2, 2, static_cast<float>(prior.variance.initial))));
    }
    else
    {
        REQUIRE(
            reader.to_map<double>("genetic/A/variance")
                .isApprox(
                    Eigen::MatrixXd::Constant(1, 2, prior.variance.initial)));
    }
    if constexpr (requires { prior.probability; })
    {
        constexpr bool sampled = requires { prior.probability.prior; };
        REQUIRE(reader.contains("genetic/A/probability") == sampled);
        if constexpr (sampled)
        {
            REQUIRE(reader.to_map<double>("genetic/A/probability")
                        .isApprox(
                            Eigen::MatrixXd::Constant(
                                1, 2, prior.probability.initial)));
        }
    }
    if constexpr (requires { prior.probabilities; })
    {
        constexpr bool sampled = requires { prior.probabilities.prior; };
        REQUIRE(reader.contains("genetic/A/probabilities") == sampled);
        if constexpr (sampled)
        {
            REQUIRE(reader.to_map<double>("genetic/A/probabilities")
                        .isApprox(
                            Eigen::Map<const Eigen::VectorXd>{
                                state.probabilities().data(),
                                gelex::ScaledMixtureState::class_count}
                                .replicate(1, 2)));
        }
        REQUIRE(reader.to_map<double>("genetic/A/component_explained_variance")
                    .isApprox(
                        Eigen::MatrixXd::Zero(
                            gelex::ScaledMixtureState::component_count, 2)));
    }
}

}  // namespace

TEST_CASE(
    "Gaussian draws preserve scalar and marker variance layouts",
    "[bayes][genetic][draws]")
{
    check_mode_draws(gelex::GaussianSpec<gelex::VarianceLayout::Pooled>{});
    check_mode_draws(gelex::GaussianSpec<gelex::VarianceLayout::Unpooled>{});
}

TEST_CASE(
    "Spike slab draws preserve assignments and optional probabilities",
    "[bayes][genetic][draws]")
{
    check_mode_draws(gelex::SpikeSlabSpec<gelex::VarianceLayout::Pooled>{});
    check_mode_draws(gelex::SpikeSlabSpec<gelex::VarianceLayout::Unpooled>{});
    check_mode_draws(
        gelex::SpikeSlabSpec<
            gelex::VarianceLayout::Pooled,
            gelex::MixtureWeightUpdate::Disabled>{});
    check_mode_draws(
        gelex::SpikeSlabSpec<
            gelex::VarianceLayout::Unpooled,
            gelex::MixtureWeightUpdate::Disabled>{});
}

TEST_CASE(
    "Scaled mixture draws preserve class payloads",
    "[bayes][genetic][draws]")
{
    check_mode_draws(gelex::ScaledMixtureSpec<>{});
    check_mode_draws(
        gelex::ScaledMixtureSpec<gelex::MixtureWeightUpdate::Disabled>{});
}

TEST_CASE(
    "Joint draws separate mode coefficients from shared assignments",
    "[bayes][genetic][draws]")
{
    const auto check = []<gelex::MixtureWeightUpdate Update>()
    {
        gelex::test::FileFixture fixture;
        const auto path = (fixture.get_test_dir() / "joint.draws").string();
        const gelex::GeneticDimensions dimensions{.individual = 3, .marker = 2};
        const gelex::HalfNormalPrior dominance_prior{
            .variance = {1.0, gelex::ScaledInvChi2LogKernel{4.0, 2.0}}};
        const gelex::JointSpikeSlabPrior<Update> joint_prior{
            .probabilities = gelex::detail::make_parameter<Update>(
                std::array<double, 4>{0.25, 0.25, 0.25, 0.25},
                gelex::make_uniform_dirichlet_prior<4>())};
        auto dominance = gelex::detail::make_state(dominance_prior, dimensions);
        auto joint = gelex::detail::make_state(joint_prior, dimensions);
        dominance.transition(0, 1.5);
        dominance.annotation_coefficients() = Eigen::Vector2d{{0.25, -0.5}};
        {
            gelex::BinaryWriter writer{path};
            auto mode_draws = gelex::detail::make_draws(
                dominance_prior, writer, "genetic/D", 1, dimensions);
            auto joint_draws = gelex::detail::make_draws(
                joint_prior, writer, "genetic/joint", 1, dimensions);
            mode_draws.append(dominance);
            joint_draws.append(joint);
            writer.close();
        }
        const gelex::BinaryReader reader{path};
        REQUIRE(reader.to_map<float>("genetic/D/coefficients")
                    .isApprox(Eigen::VectorXf{{1.5F, 0.0F}}));
        REQUIRE(reader.to_map<float>("genetic/D/annotation_coefficients")
                    .isApprox(Eigen::Vector2f{{0.25F, -0.5F}}));
        REQUIRE(reader.to_map<std::uint8_t>("genetic/joint/assignment")
                    .cast<double>()
                    .isApprox(Eigen::VectorXd::Zero(2)));
        REQUIRE(
            reader.to_map<double>("genetic/joint/component_explained_variance")
                .isApprox(Eigen::VectorXd::Zero(4)));
        REQUIRE(
            reader.contains("genetic/joint/probabilities")
            == (Update == gelex::MixtureWeightUpdate::Enabled));
    };
    check.template operator()<gelex::MixtureWeightUpdate::Enabled>();
    check.template operator()<gelex::MixtureWeightUpdate::Disabled>();
}
