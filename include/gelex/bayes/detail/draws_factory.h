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

#ifndef GELEX_BAYES_DETAIL_DRAWS_FACTORY_H_
#define GELEX_BAYES_DETAIL_DRAWS_FACTORY_H_

#include <Eigen/Core>
#include <cstddef>
#include <cstdint>
#include <fmt/format.h>
#include <string>
#include <string_view>
#include <utility>

#include "gelex/bayes/basic_draw.h"
#include "gelex/bayes/detail/state_factory.h"
#include "gelex/bayes/genetic/draws.h"
#include "gelex/bayes/genetic/prior.h"
#include "gelex/bayes/genetic/state.h"
#include "gelex/bayes/genetic_family.h"
#include "gelex/bayes/genotype/design.h"
#include "gelex/bayes/mode_values.h"
#include "gelex/genetic_mode.h"
#include "gelex/io/binary_format.h"
#include "gelex/io/binary_writer.h"

namespace gelex::detail
{

struct GeneticDrawsDimensions
{
    Eigen::Index marker_count;
    std::uint64_t draw_count;
};

class GeneticDrawsBuilder
{
   public:
    GeneticDrawsBuilder(
        BinaryWriter& writer,
        std::string prefix,
        GeneticDrawsDimensions dimensions)
        : writer_{&writer}, prefix_{std::move(prefix)}, dimensions_{dimensions}
    {
    }

    [[nodiscard]] auto scalar(std::string_view name) -> ScalarDraw
    {
        return ScalarDraw{reserve<double>(name, 1)};
    }

    [[nodiscard]] auto vector(std::string_view name, Eigen::Index rows)
        -> VectorDraw
    {
        return VectorDraw{reserve<float>(name, rows)};
    }

    template <std::size_t CategoryCount>
    [[nodiscard]] auto category(std::string_view name, Eigen::Index rows)
        -> CategoryDraw<CategoryCount>
    {
        return CategoryDraw<CategoryCount>{reserve<std::uint8_t>(name, rows)};
    }

    [[nodiscard]] auto marker_count() const noexcept -> Eigen::Index
    {
        return dimensions_.marker_count;
    }

   private:
    template <SupportedDtype T>
    [[nodiscard]] auto reserve(std::string_view name, Eigen::Index rows)
        -> PayloadWriter<T>
    {
        return writer_->reserve<T>(
            fmt::format("{}/{}", prefix_, name),
            BinaryShape{
                static_cast<std::uint64_t>(rows), dimensions_.draw_count});
    }

    BinaryWriter* writer_;
    std::string prefix_;
    GeneticDrawsDimensions dimensions_;
};

template <VarianceLayout Kind>
[[nodiscard]] auto make_marker_variance_draw(GeneticDrawsBuilder& builder)
    -> marker_variance_draw_t<Kind>
{
    if constexpr (Kind == VarianceLayout::Pooled)
    {
        return builder.scalar("variance");
    }
    else
    {
        return builder.vector("variance", builder.marker_count());
    }
}

template <UpdatePolicy Policy>
[[nodiscard]] auto make_probability_draw(
    GeneticDrawsBuilder& builder,
    std::string_view name) -> policy_draw_t<Policy, ScalarDraw>
{
    if constexpr (Policy == UpdatePolicy::Sampled)
    {
        return builder.scalar(name);
    }
    else
    {
        return EmptyDraw{};
    }
}

template <UpdatePolicy Policy, std::size_t ClassCount>
[[nodiscard]] auto make_probabilities_draw(GeneticDrawsBuilder& builder)
    -> policy_draw_t<Policy, VectorDraw>
{
    if constexpr (Policy == UpdatePolicy::Sampled)
    {
        return builder.vector(
            "probabilities", static_cast<Eigen::Index>(ClassCount));
    }
    else
    {
        return EmptyDraw{};
    }
}
// rows non-additive: they are shares, not a decomposition of genetic variance.
template <std::size_t ComponentCount>
[[nodiscard]] auto make_component_explained_variance_draw(
    GeneticDrawsBuilder& builder) -> VectorDraw
{
    return builder.vector(
        "component_explained_variance",
        static_cast<Eigen::Index>(ComponentCount));
}

template <VarianceLayout Kind>
[[nodiscard]] auto make_family_draws(
    const GaussianPrior<Kind>& /*prior*/,
    GeneticDrawsBuilder& builder) -> GaussianDraws<Kind>
{
    return {.variance = make_marker_variance_draw<Kind>(builder)};
}

[[nodiscard]] inline auto make_family_draws(
    const HalfNormalPrior<HalfNormalAsymmetry::Count>& /*prior*/,
    GeneticDrawsBuilder& builder) -> HalfNormalDraws<HalfNormalAsymmetry::Count>
{
    return {
        .variance = builder.scalar("variance"),
        .assignment = builder.category<3>("assignment", builder.marker_count()),
        .positive_probability = builder.scalar("positive_probability")};
}

[[nodiscard]] inline auto make_family_draws(
    const HalfNormalPrior<HalfNormalAsymmetry::Magnitude>& /*prior*/,
    GeneticDrawsBuilder& builder)
    -> HalfNormalDraws<HalfNormalAsymmetry::Magnitude>
{
    return {
        .variances = builder.vector("variance", 2),
        .assignment
        = builder.category<3>("assignment", builder.marker_count())};
}

template <VarianceLayout Kind, UpdatePolicy ProbabilityUpdate>
[[nodiscard]] auto make_family_draws(
    const SpikeSlabPrior<Kind, ProbabilityUpdate>& /*prior*/,
    GeneticDrawsBuilder& builder) -> SpikeSlabDraws<Kind, ProbabilityUpdate>
{
    return {
        .variance = make_marker_variance_draw<Kind>(builder),
        .assignment = builder.category<2>("assignment", builder.marker_count()),
        .probability
        = make_probability_draw<ProbabilityUpdate>(builder, "probability")};
}

template <std::size_t ClassCount, UpdatePolicy ProbabilitiesUpdate>
[[nodiscard]] auto make_family_draws(
    const ScaledMixturePrior<ClassCount, ProbabilitiesUpdate>& /*prior*/,
    GeneticDrawsBuilder& builder)
    -> ScaledMixtureDraws<ClassCount, ProbabilitiesUpdate>
{
    return {
        .variance = builder.scalar("variance"),
        .assignment
        = builder.category<ClassCount>("assignment", builder.marker_count()),
        .probabilities
        = make_probabilities_draw<ProbabilitiesUpdate, ClassCount>(builder),
        .component_explained_variance = make_component_explained_variance_draw<
            ScaledMixtureState<ClassCount>::component_count>(builder)};
}

// Rows follow the fitted column layout of JointSpikeSlabState: A in A-only,
// A in AD, D in D-only, D in AD.
template <std::size_t ClassCount, UpdatePolicy ProbabilitiesUpdate>
[[nodiscard]] auto make_family_draws(
    const JointSpikeSlabPrior<ClassCount, ProbabilitiesUpdate>& /*prior*/,
    GeneticDrawsBuilder& builder)
    -> JointSpikeSlabDraws<ClassCount, ProbabilitiesUpdate>
{
    return {
        .assignment
        = builder.category<ClassCount>("assignment", builder.marker_count()),
        .probabilities
        = make_probabilities_draw<ProbabilitiesUpdate, ClassCount>(builder),
        .component_explained_variance = make_component_explained_variance_draw<
            JointSpikeSlabState<ClassCount>::component_count>(builder)};
}

template <typename Prior>
[[nodiscard]] auto make_mode_draws(
    const Prior& prior,
    BinaryWriter& writer,
    std::string prefix,
    GeneticDrawsDimensions dimensions)
{
    GeneticDrawsBuilder builder{writer, std::move(prefix), dimensions};
    auto coefficients = builder.vector("coefficients", dimensions.marker_count);
    auto family_draws = make_family_draws(prior, builder);
    return GeneticModeDraws<decltype(family_draws)>{
        .coefficients = std::move(coefficients),
        .family_draws = std::move(family_draws)};
}

template <GeneticModeSet Modes, typename... Priors>
[[nodiscard]] auto make_genetic_draws(
    const ModeValues<Modes, Priors...>& prior,
    const bayes::GeneticDesign& design,
    BinaryWriter& writer,
    std::uint64_t draw_count)
{
    const auto dimensions = GeneticDrawsDimensions{
        .marker_count = design.cols(), .draw_count = draw_count};
    auto mode_draws = transform_mode_values(
        prior,
        [&]<GeneticMode Mode>(const auto& mode_prior)
        {
            return make_mode_draws(
                mode_prior,
                writer,
                fmt::format("genetic/{}", Mode),
                dimensions);
        });
    using GeneticState = genetic_state_t<ModeValues<Modes, Priors...>>;
    return IndependentDraws<GeneticState, decltype(mode_draws)>{
        std::move(mode_draws)};
}

template <typename ModeValuesType, typename JointPrior>
[[nodiscard]] auto make_genetic_draws(
    const JointModeValues<ModeValuesType, JointPrior>& prior,
    const bayes::GeneticDesign& design,
    BinaryWriter& writer,
    std::uint64_t draw_count)
{
    const auto dimensions = GeneticDrawsDimensions{
        .marker_count = design.cols(), .draw_count = draw_count};
    auto mode_draws = transform_mode_values(
        prior.mode_values(),
        [&]<GeneticMode Mode>(const auto& mode_prior)
        {
            return make_mode_draws(
                mode_prior,
                writer,
                fmt::format("genetic/{}", Mode),
                dimensions);
        });
    GeneticDrawsBuilder joint_builder{writer, "genetic/joint", dimensions};
    auto joint_draws = make_family_draws(prior.joint(), joint_builder);
    using GeneticState
        = genetic_state_t<JointModeValues<ModeValuesType, JointPrior>>;
    return JointDraws<
        GeneticState,
        decltype(mode_draws),
        decltype(joint_draws)>{std::move(mode_draws), std::move(joint_draws)};
}

template <typename Prior>
using genetic_draws_t = decltype(make_genetic_draws(
    std::declval<const Prior&>(),
    std::declval<const bayes::GeneticDesign&>(),
    std::declval<BinaryWriter&>(),
    std::declval<std::uint64_t>()));

}  // namespace gelex::detail

#endif  // GELEX_BAYES_DETAIL_DRAWS_FACTORY_H_
