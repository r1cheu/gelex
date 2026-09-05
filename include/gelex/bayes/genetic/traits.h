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

#ifndef GELEX_BAYES_GENETIC_TRAITS_H_
#define GELEX_BAYES_GENETIC_TRAITS_H_

#include <Eigen/Core>
#include <type_traits>

#include "gelex/bayes/basic_draw.h"
#include "gelex/bayes/basic_result.h"
#include "gelex/bayes/genetic_family.h"

namespace gelex::detail
{

template <VarianceLayout Kind>
using marker_variance_state_t = std::
    conditional_t<Kind == VarianceLayout::Pooled, double, Eigen::VectorXd>;

template <VarianceLayout Kind>
using marker_variance_draw_t = std::
    conditional_t<Kind == VarianceLayout::Pooled, ScalarDraw, VectorDraw>;

template <VarianceLayout Kind>
using marker_variance_result_t = std::conditional_t<
    Kind == VarianceLayout::Pooled,
    ScalarPosteriorResult,
    EmptyPosteriorResult>;

template <MixtureWeightUpdate Update, typename Draw>
using weight_draw_t = std::
    conditional_t<Update == MixtureWeightUpdate::Enabled, Draw, EmptyDraw>;

template <MixtureWeightUpdate Update, typename Result>
using weight_result_t = std::conditional_t<
    Update == MixtureWeightUpdate::Enabled,
    Result,
    EmptyPosteriorResult>;

}  // namespace gelex::detail

#endif  // GELEX_BAYES_GENETIC_TRAITS_H_
