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

#ifndef GELEX_BAYES_DETAIL_RESULT_FACTORY_H_
#define GELEX_BAYES_DETAIL_RESULT_FACTORY_H_

#include <Eigen/Core>
#include <cstddef>
#include <span>
#include <string>
#include <utility>
#include <vector>

#include "gelex/bayes/basic_draw.h"
#include "gelex/bayes/basic_result.h"
#include "gelex/bayes/detail/pip_factory.h"
#include "gelex/bayes/genetic/draws.h"
#include "gelex/bayes/genetic/result.h"
#include "gelex/bayes/genetic_family.h"
#include "gelex/bayes/mode_values.h"
#include "gelex/bayes/model.h"
#include "gelex/genetic_mode.h"

namespace gelex::detail
{

inline auto make_result(const EmptyDraw& /*draw*/) -> EmptyPosteriorResult
{
    return {};
}

inline auto make_result(const ScalarDraw& draw) -> ScalarPosteriorResult
{
    return ScalarPosteriorResult{std::string{draw.identifier()}, draw.result()};
}

inline auto make_result(const VectorDraw& draw) -> VectorPosteriorResult
{
    return VectorPosteriorResult{std::string{draw.identifier()}, draw.result()};
}

inline auto make_result(
    const VectorDraw& draw,
    std::span<const std::string> column_names) -> CoefficientPosteriorResult
{
    return CoefficientPosteriorResult{
        make_result(draw),
        std::vector<std::string>{column_names.begin(), column_names.end()}};
}

template <VarianceLayout Kind>
auto make_marker_variance_result(const marker_variance_draw_t<Kind>& draw)
    -> marker_variance_result_t<Kind>
{
    if constexpr (Kind == VarianceLayout::Pooled)
    {
        return make_result(draw);
    }
    else
    {
        static_cast<void>(draw);
        return {};
    }
}

template <VarianceLayout Kind>
auto make_result(const GaussianDraws<Kind>& draws)
    -> GaussianPosteriorResult<Kind>
{
    return {.variance = make_marker_variance_result<Kind>(draws.variance)};
}

inline auto make_result(const HalfNormalDraws& draws)
    -> HalfNormalPosteriorResult
{
    return {
        .variance = make_result(draws.variance),
        .probit_coefficients = make_result(draws.probit_coefficients)};
}

template <VarianceLayout Kind, MixtureWeightUpdate WeightUpdate>
auto make_result(const SpikeSlabDraws<Kind, WeightUpdate>& draws)
    -> SpikeSlabPosteriorResult<Kind, WeightUpdate>
{
    return {
        .variance = make_marker_variance_result<Kind>(draws.variance),
        .probability = make_result(draws.probability)};
}

template <std::size_t ClassCount, MixtureWeightUpdate WeightUpdate>
auto make_result(const ScaledMixtureDraws<ClassCount, WeightUpdate>& draws)
    -> ScaledMixturePosteriorResult<ClassCount, WeightUpdate>
{
    return {
        .variance = make_result(draws.variance),
        .probabilities = make_result(draws.probabilities),
        .component_explained_variance
        = make_result(draws.component_explained_variance)};
}

template <std::size_t ClassCount, MixtureWeightUpdate WeightUpdate>
auto make_result(const JointSpikeSlabDraws<ClassCount, WeightUpdate>& draws)
    -> JointSpikeSlabPosteriorResult<ClassCount, WeightUpdate>
{
    return {
        .probabilities = make_result(draws.probabilities),
        .component_explained_variance
        = make_result(draws.component_explained_variance)};
}

template <typename CoefficientDraws, typename ModeFamilyDraws>
auto make_genetic_parameters(
    const IndependentGeneticDraws<CoefficientDraws, ModeFamilyDraws>& draws)
{
    return generate_mode_values<CoefficientDraws::modes>(
        [&]<GeneticMode Mode>()
        { return make_result(draws.template family<Mode>()); });
}

template <typename CoefficientDraws, typename ModeFamilyDraws, typename JointT>
auto make_genetic_parameters(
    const JointGeneticDraws<CoefficientDraws, ModeFamilyDraws, JointT>& draws)
{
    auto mode_results = generate_mode_values<CoefficientDraws::modes>(
        [&]<GeneticMode Mode>()
        { return make_result(draws.template family<Mode>()); });
    auto joint_result = make_result(draws.joint_family());
    return JointModeValues{std::move(mode_results), std::move(joint_result)};
}

template <typename GeneticDraws>
auto make_marker_effects(const BayesModel& model, const GeneticDraws& draws)
{
    constexpr auto modes = GeneticDraws::modes;
    const double phenotype_variance = model.phenotype_variance();
    auto pip = make_pip(draws);
    auto mode_results = generate_mode_values<modes>(
        [&]<GeneticMode Mode>()
        {
            const auto& coefficients
                = draws.coefficients().template get<Mode>();
            const auto& projection = model.genetic().projection(Mode);
            Eigen::VectorXd pve = projection.col_var().transpose().array()
                                  * coefficients.mean_square().array()
                                  / phenotype_variance;
            return MarkerEffectResult{
                make_result(coefficients),
                MarkerPveResult{std::move(pve)},
                std::move(pip.template get<Mode>())};
        });

    if constexpr (modes == (GeneticMode::A | GeneticMode::D))
    {
        Eigen::VectorXd joint_pve
            = mode_results.template get<GeneticMode::A>().pve().values()
              + mode_results.template get<GeneticMode::D>().pve().values();
        const auto covariance
            = model.genetic()
                  .projection(GeneticMode::A)
                  .col_covariance(model.genetic().projection(GeneticMode::D));
        joint_pve.array() += 2.0 * covariance.transpose().array()
                             * draws.coefficients().mean_product().array()
                             / phenotype_variance;

        auto joint_pip = [&]
        {
            if constexpr (requires { pip.joint(); })
            {
                return std::move(pip.joint());
            }
            else
            {
                return EmptyPosteriorResult{};
            }
        }();
        return JointModeValues{
            std::move(mode_results),
            JointMarkerEffectResult{
                MarkerPveResult{std::move(joint_pve)}, std::move(joint_pip)}};
    }
    else
    {
        return mode_results;
    }
}

}  // namespace gelex::detail

#endif  // GELEX_BAYES_DETAIL_RESULT_FACTORY_H_
