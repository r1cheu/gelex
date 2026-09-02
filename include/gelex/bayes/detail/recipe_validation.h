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

#ifndef GELEX_BAYES_DETAIL_RECIPE_VALIDATION_H_
#define GELEX_BAYES_DETAIL_RECIPE_VALIDATION_H_

#include "gelex/bayes/variance_budget.h"
#include "gelex/types/genetic_mode.h"

namespace gelex::detail
{

auto validate_recipe_inputs(
    const VarianceBudget& variance,
    GeneticModeSet modes) -> void;

}  // namespace gelex::detail

#endif  // GELEX_BAYES_DETAIL_RECIPE_VALIDATION_H_
