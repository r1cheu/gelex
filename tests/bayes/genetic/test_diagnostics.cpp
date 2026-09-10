// Copyright 2026 RuLei Chen
// SPDX-License-Identifier: Apache-2.0

#include <Eigen/Core>
#include <array>
#include <catch2/catch_test_macros.hpp>
#include <cstddef>
#include <type_traits>
#include <variant>

#include "gelex/bayes/genetic/diagnostics_traits.h"
#include "gelex/bayes/genetic/draw_traits.h"
#include "gelex/bayes/genetic/gaussian.h"
#include "gelex/bayes/genetic/joint_spike_slab.h"
#include "gelex/bayes/genetic/parameter.h"
#include "gelex/bayes/genetic/scaled_mixture.h"
#include "gelex/bayes/genetic/spike_slab.h"
#include "gelex/bayes/genetic/types.h"
#include "gelex/bayes/mode_values.h"
#include "gelex/bayes/prior.h"
#include "gelex/bayes/recipe.h"
#include "gelex/bayes/serialization_ids.h"
#include "gelex/bayes/spec.h"
#include "gelex/bayes/stats/diagnostics.h"
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

template <typename Spec>
auto check_mode_diagnostics(const Spec& spec) -> void
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
            const double value = 0.5 * static_cast<double>(draw + 1);
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
            if constexpr (requires { state.variance() = 1.0; })
            {
                state.variance() = value;
            }
            draws << state;
        }
        sparse.close();
        dense.close();
    }
    const gelex::DenseReader dense{path};
    const gelex::CscReader sparse{sparse_path};
    const auto result = gelex::make_diagnostics(
        std::type_identity<draws_type>{},
        gelex::DrawReaders{.dense = dense, .sparse = sparse},
        gelex::genetic_id<gelex::GeneticMode::A>,
        0.9);

    constexpr bool pooled = !requires { state.variance().size(); };
    using variance_type = std::remove_cvref_t<decltype(result.variance)>;
    STATIC_REQUIRE(
        std::is_same_v<variance_type, gelex::ChainDiagnostics> == pooled);
    if constexpr (pooled)
    {
        const auto expected = gelex::diagnose_chain(
            dense.to_map<double>("genetic/A/variance").row(0), 0.9);
        REQUIRE(result.variance.mean == expected.mean);
        REQUIRE(result.variance.hpdi_upper == expected.hpdi_upper);
    }
    if constexpr (requires { prior.probability; })
    {
        constexpr bool sampled = requires { prior.probability.prior; };
        using probability_type
            = std::remove_cvref_t<decltype(result.probability)>;
        STATIC_REQUIRE(
            std::is_same_v<probability_type, gelex::ChainDiagnostics>
            == sampled);
        if constexpr (sampled)
        {
            REQUIRE(result.probability.mean == prior.probability.initial);
        }
    }
    if constexpr (requires { prior.probabilities; })
    {
        constexpr bool sampled = requires { prior.probabilities.prior; };
        using probabilities_type
            = std::remove_cvref_t<decltype(result.probabilities)>;
        STATIC_REQUIRE(
            std::is_same_v<probabilities_type, std::monostate> != sampled);
        if constexpr (sampled)
        {
            REQUIRE(
                result.probabilities.size()
                == gelex::ScaledMixtureState<>::class_count);
            REQUIRE(result.probabilities[1].mean == state.probabilities()[1]);
        }
    }
}

}  // namespace

TEST_CASE(
    "Gaussian diagnostics keep only the pooled variance",
    "[bayes][genetic][diagnostics]")
{
    check_mode_diagnostics(
        gelex::GaussianSpec<gelex::VarianceLayout::Pooled>{});
    check_mode_diagnostics(
        gelex::GaussianSpec<gelex::VarianceLayout::Unpooled>{});
}

TEST_CASE(
    "Spike slab diagnostics follow the variance and probability axes",
    "[bayes][genetic][diagnostics]")
{
    check_mode_diagnostics(
        gelex::SpikeSlabSpec<gelex::VarianceLayout::Pooled>{});
    check_mode_diagnostics(
        gelex::SpikeSlabSpec<gelex::VarianceLayout::Unpooled>{});
    check_mode_diagnostics(
        gelex::SpikeSlabSpec<
            gelex::VarianceLayout::Pooled,
            gelex::MixtureWeightUpdate::Disabled>{});
}

TEST_CASE(
    "Scaled mixture diagnostics cover variance and class probabilities",
    "[bayes][genetic][diagnostics]")
{
    check_mode_diagnostics(gelex::ScaledMixtureSpec<>{});
    check_mode_diagnostics(
        gelex::ScaledMixtureSpec<gelex::MixtureWeightUpdate::Disabled>{});
}

TEST_CASE(
    "Joint diagnostics cover half-normal and shared assignment payloads",
    "[bayes][genetic][diagnostics]")
{
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
            const double value = static_cast<double>(draw + 1);
            dominance.transition(0, value);
            dominance.variance() = value;
            dominance.annotation_coefficients()
                = Eigen::Vector2d{{value, -value}};
            joint.transition(0, 3);
            mode_draws << dominance;
            joint_draws << joint;
        }
        sparse.close();
        dense.close();
    }
    const gelex::DenseReader dense{path};
    const gelex::CscReader sparse{sparse_path};
    const gelex::DrawReaders readers{.dense = dense, .sparse = sparse};

    const auto dominance_result = gelex::make_diagnostics(
        std::type_identity<dominance_draws_type>{},
        readers,
        gelex::genetic_id<gelex::GeneticMode::D>,
        0.9);
    REQUIRE(dominance_result.variance.mean == 2.5);
    REQUIRE(dominance_result.annotation_coefficients.size() == 2);
    REQUIRE(dominance_result.annotation_coefficients[0].mean == 2.5);
    REQUIRE(dominance_result.annotation_coefficients[1].mean == -2.5);

    const auto joint_result = gelex::make_diagnostics(
        std::type_identity<joint_draws_type>{},
        readers,
        gelex::joint_genetic_id,
        0.9);
    REQUIRE(joint_result.probabilities.size() == 4);
    REQUIRE(joint_result.probabilities[3].mean == 0.25);
}
