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

#ifndef GELEX_BAYES_DETAIL_GENETIC_SPEC_H_
#define GELEX_BAYES_DETAIL_GENETIC_SPEC_H_

#include "gelex/genetic_mode.h"

namespace gelex::detail
{

template <GeneticModeSet Modes, typename Family>
struct GeneticSpecFor;

template <GeneticModeSet Modes, typename Family>
using genetic_spec_t = typename GeneticSpecFor<Modes, Family>::type;

template <GeneticModeSet Modes, typename GeneticFamily>
concept SupportedGeneticFamily
    = requires { typename genetic_spec_t<Modes, GeneticFamily>; };

}  // namespace gelex::detail

#endif  // GELEX_BAYES_DETAIL_GENETIC_SPEC_H_
