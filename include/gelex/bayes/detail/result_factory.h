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
#include <utility>

#include "gelex/bayes/basic_result.h"
#include "gelex/bayes/detail/pip_factory.h"
#include "gelex/bayes/genetic/detail/result_support.h"
#include "gelex/bayes/genetic/draws.h"
#include "gelex/bayes/genetic/gaussian.h"
#include "gelex/bayes/genetic/joint_spike_slab.h"
#include "gelex/bayes/genetic/result.h"
#include "gelex/bayes/genetic/scaled_mixture.h"
#include "gelex/bayes/genetic/spike_slab.h"
#include "gelex/bayes/mode_values.h"
#include "gelex/bayes/model.h"
#include "gelex/genetic_mode.h"

namespace gelex::detail
{

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
