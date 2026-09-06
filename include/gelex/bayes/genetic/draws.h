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

#ifndef GELEX_BAYES_GENETIC_DRAWS_H_
#define GELEX_BAYES_GENETIC_DRAWS_H_

#include <Eigen/Core>
#include <cassert>
#include <cstddef>
#include <cstdint>
#include <fmt/format.h>
#include <string>
#include <string_view>
#include <type_traits>
#include <utility>

#include "gelex/bayes/basic_draw.h"
#include "gelex/bayes/genetic/policy.h"
#include "gelex/bayes/mode_values.h"
#include "gelex/genetic_mode.h"
#include "gelex/io/binary_format.h"
#include "gelex/io/binary_writer.h"
#include "gelex/namespace.h"

GELEX_NAMESPACE_BEGIN(gelex)

template <GeneticModeSet Modes>
class GeneticCoefficientDraws
{
   public:
    static constexpr GeneticModeSet modes = Modes;

    explicit GeneticCoefficientDraws(
        HomogeneousModeValues<Modes, VectorDraw> draws)
        : draws_{std::move(draws)}
    {
    }

    template <typename GeneticState>
    auto append(const GeneticState& state) -> void
    {
        draws_.for_each(
            [&]<GeneticMode Mode>(auto& draw)
            { draw.append(state.template get<Mode>().coefficients()); });
    }

    template <GeneticMode Mode>
    [[nodiscard]] auto get() const noexcept -> const VectorDraw&
    {
        return draws_.template get<Mode>();
    }

   private:
    HomogeneousModeValues<Modes, VectorDraw> draws_;
};

template <>
class GeneticCoefficientDraws<(GeneticMode::A | GeneticMode::D)>
{
   public:
    static constexpr GeneticModeSet modes = GeneticMode::A | GeneticMode::D;

    explicit GeneticCoefficientDraws(
        HomogeneousModeValues<modes, VectorDraw> draws)
        : draws_{std::move(draws)}
    {
    }

    template <typename GeneticState>
    auto append(const GeneticState& state) -> void
    {
        const auto& additive
            = state.template get<GeneticMode::A>().coefficients();
        const auto& dominance
            = state.template get<GeneticMode::D>().coefficients();

        draws_.template get<GeneticMode::A>().append(additive);
        draws_.template get<GeneticMode::D>().append(dominance);

        const Eigen::ArrayXd product = additive.array() * dominance.array();
        ++count_;
        if (count_ == 1)
        {
            mean_product_ = product;
        }
        else
        {
            mean_product_.array() += (product - mean_product_.array())
                                     / static_cast<double>(count_);
        }
    }

    template <GeneticMode Mode>
    [[nodiscard]] auto get() const noexcept -> const VectorDraw&
    {
        return draws_.template get<Mode>();
    }

    [[nodiscard]] auto mean_product() const noexcept -> const Eigen::VectorXd&
    {
        assert(count_ > 0);
        return mean_product_;
    }

   private:
    HomogeneousModeValues<modes, VectorDraw> draws_;
    Eigen::VectorXd mean_product_;
    std::uint64_t count_{0};
};

template <typename CoefficientDraws, typename ModeFamilyDraws>
class IndependentGeneticDraws
{
   public:
    static constexpr GeneticModeSet modes = CoefficientDraws::modes;

    IndependentGeneticDraws(
        CoefficientDraws coefficients,
        ModeFamilyDraws families)
        : coefficients_{std::move(coefficients)}, families_{std::move(families)}
    {
    }

    template <typename GeneticState>
    auto append(const GeneticState& state) -> void
    {
        coefficients_.append(state);
        families_.for_each([&]<GeneticMode Mode>(auto& family)
                           { family.append(state.template get<Mode>()); });
    }

    [[nodiscard]] auto coefficients() const noexcept -> const CoefficientDraws&
    {
        return coefficients_;
    }

    template <GeneticMode Mode>
    [[nodiscard]] auto family() const noexcept -> const
        typename ModeFamilyDraws::template mode_value_type<Mode>&
    {
        return families_.template get<Mode>();
    }

   private:
    CoefficientDraws coefficients_;
    ModeFamilyDraws families_;
};

template <typename CoefficientDraws, typename ModeFamilyDraws, typename JointT>
class JointGeneticDraws
{
   public:
    static constexpr GeneticModeSet modes = CoefficientDraws::modes;

    JointGeneticDraws(
        CoefficientDraws coefficients,
        ModeFamilyDraws families,
        JointT joint_family)
        : coefficients_{std::move(coefficients)},
          families_{std::move(families)},
          joint_family_{std::move(joint_family)}
    {
    }

    template <typename GeneticState>
    auto append(const GeneticState& state) -> void
    {
        coefficients_.append(state);
        families_.for_each([&]<GeneticMode Mode>(auto& family)
                           { family.append(state.template get<Mode>()); });
        joint_family_.append(state.joint());
    }

    [[nodiscard]] auto coefficients() const noexcept -> const CoefficientDraws&
    {
        return coefficients_;
    }

    template <GeneticMode Mode>
    [[nodiscard]] auto family() const noexcept -> const
        typename ModeFamilyDraws::template mode_value_type<Mode>&
    {
        return families_.template get<Mode>();
    }

    [[nodiscard]] auto joint_family() const noexcept -> const JointT&
    {
        return joint_family_;
    }

   private:
    CoefficientDraws coefficients_;
    ModeFamilyDraws families_;
    JointT joint_family_;
};

GELEX_NAMESPACE_BEGIN(detail)
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
using marker_variance_draw_t = std::
    conditional_t<Kind == VarianceLayout::Pooled, ScalarDraw, VectorDraw>;

template <MixtureWeightUpdate Update, typename Draw>
using weight_draw_t = std::
    conditional_t<Update == MixtureWeightUpdate::Enabled, Draw, EmptyDraw>;

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

template <MixtureWeightUpdate Update>
[[nodiscard]] auto make_probability_draw(
    GeneticDrawsBuilder& builder,
    std::string_view name) -> weight_draw_t<Update, ScalarDraw>
{
    if constexpr (Update == MixtureWeightUpdate::Enabled)
    {
        return builder.scalar(name);
    }
    else
    {
        return EmptyDraw{};
    }
}

template <MixtureWeightUpdate Update, std::size_t ClassCount>
[[nodiscard]] auto make_probabilities_draw(GeneticDrawsBuilder& builder)
    -> weight_draw_t<Update, VectorDraw>
{
    if constexpr (Update == MixtureWeightUpdate::Enabled)
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

[[nodiscard]] constexpr auto is_non_null_category(std::size_t category) noexcept
    -> bool
{
    return category != 0;
}
GELEX_NAMESPACE_END(detail)

GELEX_NAMESPACE_END(gelex)

#endif  // GELEX_BAYES_GENETIC_DRAWS_H_
