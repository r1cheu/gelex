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

#ifndef GELEX_BAYES_GENETIC_DETAIL_STATE_SUPPORT_H_
#define GELEX_BAYES_GENETIC_DETAIL_STATE_SUPPORT_H_

#include <Eigen/Core>

#include "gelex/bayes/genetic/traits.h"
#include "gelex/bayes/genetic_family.h"
#include "gelex/bayes/parameter.h"

namespace gelex::detail
{

struct GeneticStateDimensions
{
    Eigen::Index marker_count;
    Eigen::Index individual_count;
};

template <VarianceLayout Kind>
auto initial_marker_variance(
    const VarianceParameter& parameter,
    Eigen::Index marker_count) -> marker_variance_state_t<Kind>
{
    if constexpr (Kind == VarianceLayout::Pooled)
    {
        return parameter.initial;
    }
    else
    {
        return Eigen::VectorXd::Constant(marker_count, parameter.initial);
    }
}

}  // namespace gelex::detail

#endif  // GELEX_BAYES_GENETIC_DETAIL_STATE_SUPPORT_H_
