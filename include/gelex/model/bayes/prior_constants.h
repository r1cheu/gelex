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

#ifndef GELEX_MODEL_BAYES_PRIOR_CONSTANTS_H_
#define GELEX_MODEL_BAYES_PRIOR_CONSTANTS_H_

namespace gelex
{

// Prior distribution constants for Bayesian models
namespace prior_constants
{

// Inverse-Gamma prior parameters for marker variance
// Shape parameter for marker variance prior
constexpr double MARKER_VARIANCE_SHAPE = 4.0;
// Scale parameter multiplier for marker variance prior
constexpr double MARKER_VARIANCE_SCALE_MULTIPLIER = 0.5;

// Inverse-Gamma prior parameters for residual variance
// Shape parameter for residual variance prior
constexpr double RESIDUAL_VARIANCE_SHAPE = 4.0;
// Scale parameter for residual variance prior
constexpr double RESIDUAL_VARIANCE_SCALE = 0.0;

// Inverse-Gamma prior parameters for random effects variance
// Shape parameter for random effects variance prior
constexpr double RANDOM_EFFECTS_SHAPE = 4.0;
// Scale parameter for random effects variance prior
constexpr double RANDOM_EFFECTS_SCALE = 0.0;

// Default non-zero marker proportion for non mixture model
constexpr double NON_MIXTURE_PROPORTION = 1.0;
}  // namespace prior_constants

}  // namespace gelex

#endif  // GELEX_MODEL_BAYES_PRIOR_CONSTANTS_H_
