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
#include <utility>

#include "gelex/bayes/basic_draw.h"
#include "gelex/bayes/genetic/state.h"
#include "gelex/bayes/genetic/traits.h"
#include "gelex/bayes/genetic_family.h"
#include "gelex/bayes/mode_values.h"
#include "gelex/genetic_mode.h"
#include "gelex/infra/var.h"

namespace gelex
{

struct HalfNormalDraws
{
    ScalarDraw variance;
    VectorDraw probit_coefficients;

    auto append(const HalfNormalState& state) -> void
    {
        variance.append(state.variance);
        probit_coefficients.append(state.probit_coefficients);
    }
};

template <VarianceLayout Kind, MixtureWeightUpdate WeightUpdate>
struct SpikeSlabDraws
{
    detail::marker_variance_draw_t<Kind> variance;
    CategoryDraw<2> assignment;
    detail::weight_draw_t<WeightUpdate, ScalarDraw> probability;

    auto append(const SpikeSlabState<Kind>& state) -> void
    {
        variance.append(state.variance);
        assignment.append(state.assignment);
        probability.append(state.probability);
    }
};

template <std::size_t ClassCount, MixtureWeightUpdate WeightUpdate>
struct ScaledMixtureDraws
{
    ScalarDraw variance;
    CategoryDraw<ClassCount> assignment;
    detail::weight_draw_t<WeightUpdate, VectorDraw> probabilities;
    VectorDraw component_explained_variance;

    auto append(const ScaledMixtureState<ClassCount>& state) -> void
    {
        variance.append(state.variance);
        assignment.append(state.assignment);
        probabilities.append(state.probabilities);
        component_explained_variance.append(
            matvar<0>(state.fitted_values, VarNormType::Population));
    }
};

template <std::size_t ClassCount, MixtureWeightUpdate WeightUpdate>
struct JointSpikeSlabDraws
{
    CategoryDraw<ClassCount> assignment;
    detail::weight_draw_t<WeightUpdate, VectorDraw> probabilities;
    VectorDraw component_explained_variance;

    auto append(const JointSpikeSlabState<ClassCount>& state) -> void
    {
        assignment.append(state.assignment);
        probabilities.append(state.probabilities);
        component_explained_variance.append(
            matvar<0>(state.fitted_values, VarNormType::Population));
    }
};

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
            { draw.append(state.template get<Mode>().coefficients); });
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
            = state.template get<GeneticMode::A>().coefficients;
        const auto& dominance
            = state.template get<GeneticMode::D>().coefficients;

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
        families_.for_each(
            [&]<GeneticMode Mode>(auto& family)
            { family.append(state.template get<Mode>().family_state); });
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
        families_.for_each(
            [&]<GeneticMode Mode>(auto& family)
            {
                family.append(
                    state.mode_values().template get<Mode>().family_state);
            });
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

}  // namespace gelex

#endif  // GELEX_BAYES_GENETIC_DRAWS_H_
