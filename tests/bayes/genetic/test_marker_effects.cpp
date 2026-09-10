// Copyright 2026 RuLei Chen
// SPDX-License-Identifier: Apache-2.0

#include <Eigen/Core>
#include <array>
#include <catch2/catch_test_macros.hpp>
#include <cstddef>
#include <string>
#include <type_traits>
#include <utility>
#include <vector>

#include "gelex/bayes/genetic/diagnostics_traits.h"
#include "gelex/bayes/genetic/draw_traits.h"
#include "gelex/bayes/genetic/gaussian.h"
#include "gelex/bayes/genetic/joint_spike_slab.h"
#include "gelex/bayes/genetic/marker_effect_traits.h"
#include "gelex/bayes/genetic/parameter.h"
#include "gelex/bayes/genetic/scaled_mixture.h"
#include "gelex/bayes/genetic/spike_slab.h"
#include "gelex/bayes/genetic/types.h"
#include "gelex/bayes/mode_values.h"
#include "gelex/bayes/prior.h"
#include "gelex/bayes/recipe.h"
#include "gelex/bayes/serialization_ids.h"
#include "gelex/bayes/spec.h"
#include "gelex/bayes/stats/dirichlet_log_kernel.h"
#include "gelex/bayes/stats/scaled_inv_chi2_log_kernel.h"
#include "gelex/bayes/variance/budget.h"
#include "gelex/genetic_mode.h"
#include "gelex/io/csc_reader.h"
#include "gelex/io/csc_writer.h"
#include "gelex/io/dense_reader.h"
#include "gelex/io/dense_writer.h"

#include "bayes/bayes_model_fixture.h"
#include "file_fixture.h"

namespace
{

constexpr std::size_t draw_count = 4;

auto names_of(const gelex::MarkerEffectTable& table) -> std::vector<std::string>
{
    return {table.names().begin(), table.names().end()};
}

// Marker 0 takes draws 0.5, 1.0, 1.5, 2.0 (mixture families: always active);
// marker 1 stays at zero.
template <typename Spec>
auto check_mode_marker_effects(const Spec& spec) -> void
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
    const auto sparse_path = path + ".csc";
    const gelex::GeneticDimensions dimensions{.individual = 3, .marker = 2};
    auto state = gelex::make_state(prior, dimensions);
    using draws_type = decltype(gelex::make_draws(
        state,
        std::declval<gelex::DrawWriters>(),
        gelex::genetic_id<gelex::GeneticMode::A>,
        draw_count));
    {
        gelex::DenseWriter dense{path};
        gelex::CscWriter sparse{sparse_path};
        auto draws = gelex::make_draws(
            state,
            gelex::DrawWriters{.dense = dense, .sparse = sparse},
            gelex::genetic_id<gelex::GeneticMode::A>,
            draw_count);
        for (std::size_t draw = 0; draw < draw_count; ++draw)
        {
            const auto value = 0.5 * static_cast<double>(draw + 1);
            if constexpr (requires { state.probability(); })
            {
                state.transition(0, value, true);
            }
            else if constexpr (requires { state.assignments(); })
            {
                state.transition(0, value, 1);
            }
            else
            {
                state.transition(0, value);
            }
            draws << state;
        }
        sparse.close();
        dense.close();
    }
    const gelex::DenseReader dense{path};
    const gelex::CscReader sparse{sparse_path};
    const auto result = gelex::make_marker_effects(
        std::type_identity<draws_type>{},
        gelex::DrawReaders{.dense = dense, .sparse = sparse},
        gelex::genetic_id<gelex::GeneticMode::A>);

    REQUIRE(result.coefficients.mean.isApprox(Eigen::VectorXd{{1.25, 0.0}}));
    REQUIRE(result.coefficients.sd(0) > 0.0);
    REQUIRE(result.coefficients.sd(1) == 0.0);

    gelex::MarkerEffectTable table{2};
    gelex::append_marker_columns(
        table,
        result,
        gelex::GeneticMode::A,
        gelex::MarkerPveScale{Eigen::RowVectorXd{{1.0, 1.0}}, 1.0});
    std::vector<std::string> expected{"BETA_A", "SE_A", "PVE_A"};
    if constexpr (requires { result.pip; })
    {
        REQUIRE(result.pip.isApprox(Eigen::VectorXd{{1.0, 0.0}}));
        expected.emplace_back("PIP_A");
    }
    REQUIRE(names_of(table) == expected);
    REQUIRE(table.column("BETA_A").isApprox(result.coefficients.mean));
}

}  // namespace

TEST_CASE(
    "Gaussian marker effects carry coefficient summaries only",
    "[bayes][genetic][marker_effects]")
{
    check_mode_marker_effects(
        gelex::GaussianSpec<gelex::VarianceLayout::Pooled>{});
    check_mode_marker_effects(
        gelex::GaussianSpec<gelex::VarianceLayout::Unpooled>{});
}

TEST_CASE(
    "Spike slab marker effects add the inclusion probability",
    "[bayes][genetic][marker_effects]")
{
    check_mode_marker_effects(
        gelex::SpikeSlabSpec<gelex::VarianceLayout::Pooled>{});
    check_mode_marker_effects(
        gelex::SpikeSlabSpec<
            gelex::VarianceLayout::Unpooled,
            gelex::MixtureWeightUpdate::Disabled>{});
}

TEST_CASE(
    "Scaled mixture marker effects add the inclusion probability",
    "[bayes][genetic][marker_effects]")
{
    check_mode_marker_effects(gelex::ScaledMixtureSpec<>{});
}

TEST_CASE(
    "Joint marker effects split the shared assignment by mode",
    "[bayes][genetic][marker_effects]")
{
    STATIC_REQUIRE(gelex::joint_class_activates(1, gelex::GeneticMode::A));
    STATIC_REQUIRE(!gelex::joint_class_activates(1, gelex::GeneticMode::D));
    STATIC_REQUIRE(gelex::joint_class_activates(2, gelex::GeneticMode::D));
    STATIC_REQUIRE(gelex::joint_class_activates(3, gelex::GeneticMode::A));
    STATIC_REQUIRE(gelex::joint_class_activates(3, gelex::GeneticMode::D));
    STATIC_REQUIRE(!gelex::joint_class_activates(0, gelex::GeneticMode::A));

    gelex::test::FileFixture fixture;
    const auto path = (fixture.get_test_dir() / "joint.draws").string();
    const auto sparse_path = path + ".csc";
    const gelex::GeneticDimensions dimensions{.individual = 3, .marker = 2};
    const gelex::HalfNormalPrior dominance_prior{
        .variance = {1.0, gelex::ScaledInvChi2LogKernel{4.0, 2.0}}};
    const gelex::JointSpikeSlabPrior<> joint_prior{
        .probabilities
        = gelex::detail::make_parameter<gelex::MixtureWeightUpdate::Enabled>(
            std::array<double, 4>{0.25, 0.25, 0.25, 0.25},
            gelex::make_uniform_dirichlet_prior<4>())};
    auto dominance = gelex::make_state(dominance_prior, dimensions);
    auto joint = gelex::make_state(joint_prior, dimensions);
    using dominance_draws_type
        = gelex::HalfNormalDraws<gelex::CoefficientLayout::Sparse>;
    using joint_draws_type
        = gelex::JointSpikeSlabDraws<gelex::MixtureWeightUpdate::Enabled>;
    {
        gelex::DenseWriter dense{path};
        gelex::CscWriter sparse{sparse_path};
        const gelex::DrawWriters writers{.dense = dense, .sparse = sparse};
        auto mode_draws = gelex::make_draws<gelex::CoefficientLayout::Sparse>(
            dominance,
            writers,
            gelex::genetic_id<gelex::GeneticMode::D>,
            draw_count);
        auto joint_draws = gelex::make_draws(
            joint, writers, gelex::joint_genetic_id, draw_count);
        for (std::size_t draw = 0; draw < draw_count; ++draw)
        {
            const auto value = static_cast<double>(draw + 1);
            dominance.transition(0, value);
            // marker 0: A, A, AD, AD; marker 1: D in the last draw only
            joint.transition(0, draw < 2 ? 1 : 3);
            joint.transition(1, draw == 3 ? 2 : 0);
            mode_draws << dominance;
            joint_draws << joint;
        }
        sparse.close();
        dense.close();
    }
    const gelex::DenseReader dense{path};
    const gelex::CscReader sparse{sparse_path};
    const gelex::DrawReaders readers{.dense = dense, .sparse = sparse};

    const auto dominance_result = gelex::make_marker_effects(
        std::type_identity<dominance_draws_type>{},
        readers,
        gelex::genetic_id<gelex::GeneticMode::D>);
    REQUIRE(dominance_result.coefficients.mean.isApprox(
        Eigen::VectorXd{{2.5, 0.0}}));

    const auto joint_result = gelex::make_marker_effects(
        std::type_identity<joint_draws_type>{},
        readers,
        gelex::joint_genetic_id);
    REQUIRE(joint_result.mode_pip.get<gelex::GeneticMode::A>().isApprox(
        Eigen::VectorXd{{1.0, 0.0}}));
    REQUIRE(joint_result.mode_pip.get<gelex::GeneticMode::D>().isApprox(
        Eigen::VectorXd{{0.5, 0.25}}));
    REQUIRE(joint_result.pip.isApprox(Eigen::VectorXd{{1.0, 0.25}}));

    gelex::MarkerEffectTable table{2};
    gelex::append_marker_columns(table, joint_result, gelex::GeneticMode::D);
    gelex::append_marker_columns(table, joint_result);
    REQUIRE(names_of(table) == std::vector<std::string>{"PIP_D", "PIP"});
    REQUIRE(table.column("PIP_D")(1) == 0.25);
}
