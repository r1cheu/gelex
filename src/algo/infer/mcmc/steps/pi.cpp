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

#include "gelex/algo/infer/mcmc/steps/pi.h"

#include <utility>

#include "gelex/algo/infer/mcmc/context.h"
#include "gelex/exception.h"
#include "gelex/model/bayes/prior_parameters.h"
#include "gelex/types/genetic_effect_type.h"

namespace gelex::mcmc
{

PiStep::PiStep(Deps deps) : deps_(std::move(deps)), dirichlet_(deps_.alpha) {}

auto PiStep::make(const Context&, GeneticMode) -> PiStep
{
    throw GelexException(
        "PiStep is not implemented after Bayes prior/state cleanup");
}

auto PiStep::step() -> void
{
    dirichlet_.reset();
    if (deps_.proportion.update == bayes::UpdatePolicy::sampled)
    {
        deps_.proportion.value = dirichlet_(deps_.proportion.count, deps_.rng);
    }
}

}  // namespace gelex::mcmc
