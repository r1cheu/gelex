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

#ifndef GELEX_BAYES_GENETIC_PARAMETER_H_
#define GELEX_BAYES_GENETIC_PARAMETER_H_

#include <array>
#include <cstddef>
#include <type_traits>

#include "gelex/bayes/genetic_family.h"
#include "gelex/bayes/parameter.h"
#include "gelex/bayes/stats/dirichlet_log_kernel.h"

namespace gelex
{

template <MixtureWeightUpdate Update>
using ProbabilityParameter = std::conditional_t<
    Update == MixtureWeightUpdate::Disabled,
    FixedParameter<double>,
    Parameter<double, DirichletLogKernel<2>>>;

template <std::size_t ClassCount, MixtureWeightUpdate Update>
using SimplexParameter = std::conditional_t<
    Update == MixtureWeightUpdate::Disabled,
    FixedParameter<std::array<double, ClassCount>>,
    Parameter<std::array<double, ClassCount>, DirichletLogKernel<ClassCount>>>;

}  // namespace gelex

#endif  // GELEX_BAYES_GENETIC_PARAMETER_H_
