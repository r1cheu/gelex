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

#include "gelex/algo/infer/vi/updaters/additive.h"

#include "gelex/algo/infer/vi/updaters/rr.h"
#include "gelex/model/bayes/model.h"
#include "gelex/model/bayes/prior.h"

namespace gelex::detail::AdditiveUpdater
{

auto RR::operator()(
    const BayesModel& model,
    const bayes::Priors& priors,
    vi::State& states) const -> void
{
    gelex::vi::detail::RR(
        *model.genetic(GeneticMode::A),
        *priors.genetic(GeneticMode::A),
        *states.genetic(GeneticMode::A),
        states.residual());
}

}  // namespace gelex::detail::AdditiveUpdater
