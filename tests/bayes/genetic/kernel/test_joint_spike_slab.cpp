// Copyright 2026 RuLei Chen
// SPDX-License-Identifier: Apache-2.0

#include <Eigen/Core>
#include <array>
#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>
#include <concepts>
#include <cstddef>
#include <cstdint>
#include <random>
#include <utility>
#include <vector>

#include "gelex/bayes/genetic/gaussian.h"
#include "gelex/bayes/genetic/joint_spike_slab.h"
#include "gelex/bayes/genetic/marker_covariate.h"
#include "gelex/bayes/genetic/marker_covariate_io.h"
#include "gelex/bayes/genetic/types.h"
#include "gelex/bayes/genotype/design.h"
#include "gelex/bayes/genotype/projection.h"
#include "gelex/bayes/kernel.h"
#include "gelex/bayes/mode_values.h"
#include "gelex/bayes/model.h"
#include "gelex/bayes/prior.h"
#include "gelex/bayes/recipe.h"
#include "gelex/bayes/spec.h"
#include "gelex/bayes/state.h"
#include "gelex/bayes/variance/budget.h"
#include "gelex/data/bed.h"
#include "gelex/data/fixed_design.h"
#include "gelex/data/genotype_method.h"
#include "gelex/genetic_mode.h"

#include "bed_fixture.h"

using Catch::Approx;

namespace
{

constexpr auto mode_ad = gelex::GeneticMode::A | gelex::GeneticMode::D;

using JointSpikeSlabSpecAD = gelex::JointModeValues<
    gelex::ModeValues<mode_ad, gelex::GaussianSpec<>, gelex::HalfNormalSpec>,
    gelex::JointSpikeSlabSpec<>>;
using FixedJointSpikeSlabSpecAD = gelex::JointModeValues<
    gelex::ModeValues<mode_ad, gelex::GaussianSpec<>, gelex::HalfNormalSpec>,
    gelex::JointSpikeSlabSpec<gelex::MixtureWeightUpdate::Disabled>>;

using JointModeSpecs
    = gelex::ModeValues<mode_ad, gelex::GaussianSpec<>, gelex::HalfNormalSpec>;

using JointModePriors = gelex::ModeValues<
    mode_ad,
    gelex::GaussianPrior<gelex::VarianceLayout::Pooled>,
    gelex::HalfNormalPrior>;

template <gelex::MixtureWeightUpdate Update>
using JointGeneticPrior = gelex::
    JointModeValues<JointModePriors, gelex::JointSpikeSlabPrior<Update>>;

using SampledPrior = JointGeneticPrior<gelex::MixtureWeightUpdate::Enabled>;
using FixedPrior = JointGeneticPrior<gelex::MixtureWeightUpdate::Disabled>;

struct AssignmentObservation
{
    std::size_t class_index{};
    double coefficient{};
};

static_assert(std::same_as<
              decltype(gelex::make_kernel(
                  std::declval<const gelex::BayesPrior<SampledPrior>&>())),
              gelex::BayesKernel<SampledPrior>>);
static_assert(std::same_as<
              decltype(gelex::make_kernel(
                  std::declval<const gelex::BayesPrior<FixedPrior>&>())),
              gelex::BayesKernel<FixedPrior>>);

auto make_model() -> gelex::BayesModel
{
    const Eigen::MatrixXd genotypes{
        {0.0, 0.0, 1.0, 2.0},
        {1.0, 0.0, 2.0, 1.0},
        {2.0, 1.0, 0.0, 0.0},
        {0.0, 2.0, 1.0, 1.0},
        {1.0, 2.0, 0.0, 2.0},
        {2.0, 1.0, 2.0, 0.0}};
    gelex::test::BedFixture fixture;
    const auto [prefix, raw] = fixture.create_deterministic_bed_files(
        genotypes,
        {},
        {"marker_1", "marker_2", "marker_3", "marker_4"},
        {"1", "1", "1", "1"},
        {{'A', 'G'}, {'A', 'G'}, {'A', 'G'}, {'A', 'G'}});
    static_cast<void>(raw);
    auto bed = gelex::open_bed(prefix.string());

    const auto annotation_path = fixture.get_file_fixture().create_text_file(
        "CHR\tSNP\tBP\tA1\tA2\tAnnotation\n"
        "1\tmarker_1\t1\tA\tG\t-0.5\n"
        "1\tmarker_2\t2\tA\tG\t0.0\n"
        "1\tmarker_3\t3\tA\tG\t0.5\n"
        "1\tmarker_4\t4\tA\tG\t1.0\n",
        ".anno");
    auto annotation_frame
        = gelex::bayes::read_marker_annotation(annotation_path);
    auto marker_covariate = gelex::bayes::make_marker_covariate(
        std::move(annotation_frame), bed.bim());
    auto genetic = gelex::bayes::GeneticDesign{
        std::move(bed),
        mode_ad,
        gelex::GenotypeMethod::NOIACenter,
        std::move(marker_covariate)};
    return gelex::BayesModel{
        Eigen::VectorXd{{1.0, -0.5, 0.25, 2.0, -1.0, 0.75}},
        gelex::FixedDesign::make(genotypes.rows()),
        {},
        std::move(genetic)};
}

auto reconstruct_total(
    const gelex::bayes::GeneticDesign& design,
    gelex::GeneticMode mode,
    const Eigen::VectorXd& coefficients) -> Eigen::VectorXd
{
    Eigen::VectorXd fitted_values = Eigen::VectorXd::Zero(design.rows());
    const auto& projection = design.projection(mode);
    for (Eigen::Index marker = 0; marker < coefficients.size(); ++marker)
    {
        projection.axpy(marker, coefficients(marker), fitted_values);
    }
    return fitted_values;
}

template <typename State>
auto initialize_non_null_state(const gelex::BayesModel& model, State& state)
    -> void
{
    auto& mode_states = state.genetic().mode_values();
    auto& additive = mode_states.template get<gelex::GeneticMode::A>();
    auto& dominance = mode_states.template get<gelex::GeneticMode::D>();
    auto& joint = state.genetic().joint();

    const Eigen::VectorXd additive_values{{0.0, 0.3, 0.0, -0.2}};
    const Eigen::VectorXd dominance_values{{0.0, 0.0, -0.4, 0.5}};
    for (Eigen::Index marker = 0; marker < 4; ++marker)
    {
        joint.transition(marker, static_cast<std::uint8_t>(marker));
        additive.transition(marker, additive_values(marker));
        dominance.transition(marker, dominance_values(marker));
    }
    state.residual().adjusted_response -= reconstruct_total(
        model.genetic(), gelex::GeneticMode::A, additive.coefficients());
    state.residual().adjusted_response -= reconstruct_total(
        model.genetic(), gelex::GeneticMode::D, dominance.coefficients());
}

template <typename State>
auto require_residual_invariant(
    const gelex::BayesModel& model,
    const State& state) -> void
{
    const auto& mode_states = state.genetic().mode_values();
    const auto& additive = mode_states.template get<gelex::GeneticMode::A>();
    const auto& dominance = mode_states.template get<gelex::GeneticMode::D>();
    const auto additive_fitted_values = reconstruct_total(
        model.genetic(), gelex::GeneticMode::A, additive.coefficients());
    const auto dominance_fitted_values = reconstruct_total(
        model.genetic(), gelex::GeneticMode::D, dominance.coefficients());
    const Eigen::VectorXd fixed_fitted_values
        = model.fixed().X() * state.fixed().coefficients;

    REQUIRE((state.residual().adjusted_response + fixed_fitted_values
             + additive_fitted_values + dominance_fitted_values)
                .isApprox(model.phenotype()));
}

auto require_additive_assignment(const AssignmentObservation& observation)
    -> void
{
    const bool active
        = observation.class_index == 1 || observation.class_index == 3;
    if (!active)
    {
        REQUIRE(observation.coefficient == 0.0);
    }
}

auto require_dominance_assignment(const AssignmentObservation& observation)
    -> void
{
    const bool active
        = observation.class_index == 2 || observation.class_index == 3;
    if (active)
    {
        REQUIRE(observation.coefficient != 0.0);
        return;
    }
    REQUIRE(observation.coefficient == 0.0);
}

template <typename State>
auto require_assignment_invariants(const State& state) -> void
{
    const auto& mode_states = state.genetic().mode_values();
    const auto& additive = mode_states.template get<gelex::GeneticMode::A>();
    const auto& dominance = mode_states.template get<gelex::GeneticMode::D>();
    const auto& joint = state.genetic().joint();

    for (Eigen::Index marker = 0; marker < joint.assignments().size(); ++marker)
    {
        const auto class_index
            = static_cast<std::size_t>(joint.assignments()(marker));
        REQUIRE(class_index < gelex::JointSpikeSlabSpec<>::class_count);
        require_additive_assignment(
            {.class_index = class_index,
             .coefficient = additive.coefficients()(marker)});
        require_dominance_assignment(
            {.class_index = class_index,
             .coefficient = dominance.coefficients()(marker)});
    }
}

template <typename State>
auto require_state_invariants(
    const gelex::BayesModel& model,
    const State& state) -> void
{
    require_residual_invariant(model, state);
    require_assignment_invariants(state);
}

template <typename Probabilities>
auto require_probability_simplex(const Probabilities& probabilities) -> void
{
    double sum = 0.0;
    for (const double probability : probabilities)
    {
        REQUIRE(probability > 0.0);
        sum += probability;
    }
    REQUIRE(sum == Approx(1.0));
}

TEST_CASE(
    "joint half-normal kernel maintains totals and fixed component groups",
    "[bayes][kernel][joint_spike_slab]")
{
    constexpr std::array probabilities{0.25, 0.25, 0.25, 0.25};
    const auto model = make_model();
    const auto prior = gelex::make_prior(
        gelex::BayesRecipe<mode_ad, JointSpikeSlabSpecAD>{
            JointSpikeSlabSpecAD{
                JointModeSpecs{
                    gelex::GaussianSpec<>{}, gelex::HalfNormalSpec{}},
                gelex::JointSpikeSlabSpec<>{probabilities}},
            gelex::VarianceBudget{{.additive = 0.4, .dominance = 0.1}}},
        model);
    auto state = gelex::make_state(prior, model);
    auto kernel = gelex::make_kernel(prior);
    std::mt19937_64 rng{123};
    initialize_non_null_state(model, state);

    kernel.step(model, state, rng);

    require_state_invariants(model, state);
    const auto& mode_states = state.genetic().mode_values();
    const auto& additive = mode_states.get<gelex::GeneticMode::A>();
    const auto& dominance = mode_states.get<gelex::GeneticMode::D>();
    REQUIRE(additive.variance() > 0.0);
    REQUIRE(dominance.variance() > 0.0);
    REQUIRE(dominance.annotation_coefficients().allFinite());
    require_probability_simplex(state.genetic().joint().probabilities());
}

TEST_CASE(
    "joint fixed allocation kernel preserves class probabilities",
    "[bayes][kernel][joint_spike_slab]")
{
    constexpr std::array probabilities{0.1, 0.2, 0.3, 0.4};
    const auto model = make_model();
    const auto prior = gelex::make_prior(
        gelex::BayesRecipe<mode_ad, FixedJointSpikeSlabSpecAD>{
            FixedJointSpikeSlabSpecAD{
                JointModeSpecs{
                    gelex::GaussianSpec<>{}, gelex::HalfNormalSpec{}},
                gelex::JointSpikeSlabSpec<gelex::MixtureWeightUpdate::Disabled>{
                    probabilities}},
            gelex::VarianceBudget{{.additive = 0.4, .dominance = 0.1}}},
        model);
    auto state = gelex::make_state(prior, model);
    auto kernel = gelex::make_kernel(prior);
    std::mt19937_64 rng{321};
    initialize_non_null_state(model, state);

    kernel.step(model, state, rng);

    require_state_invariants(model, state);
    REQUIRE(state.genetic().joint().probabilities() == probabilities);
}

}  // namespace
