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

#ifndef GELEX_INFRA_LOGGING_FIT_EVENT_H_
#define GELEX_INFRA_LOGGING_FIT_EVENT_H_

#include <cstddef>
#include <functional>
#include <string>
#include <variant>
#include <vector>

#include <Eigen/Core>

#include "gelex/algo/infer/mcmc/state.h"
#include "gelex/model/bayes/model.h"
#include "gelex/model/bayes/prior.h"
#include "gelex/model/bayes/recipe.h"
#include "gelex/types/genetic_effect_type.h"

namespace gelex
{

namespace mcmc
{
class Result;
}  // namespace mcmc

// --- MCMC-specific events ---

struct MCMCBannerEvent
{
};

struct MCMCConfigEvent
{
    gelex::bayes::BayesRecipePreset preset;
    std::vector<GeneticMode> requested_effects;
    Eigen::Index n_iters{};
    Eigen::Index n_burn_in{};
    int seed{};
};

struct MCMCProgressEvent
{
    size_t current{};
    size_t total{};
    bool done{};
    const mcmc::State* state{};
};

struct MCMCCompleteEvent
{
    const mcmc::Result* result;
    const BayesModel* model;
    std::ptrdiff_t samples_collected;
};

// --- Shared fit protocol events ---

struct FitPriorSetEvent
{
    const bayes::BayesPrior* prior{};
};

struct FitResultsSavedEvent
{
    std::string out_prefix;
};

struct FitCheckpointSavedEvent
{
};

using MCMCEvent = std::variant<
    MCMCBannerEvent,
    MCMCConfigEvent,
    FitPriorSetEvent,
    MCMCProgressEvent,
    MCMCCompleteEvent,
    FitResultsSavedEvent,
    FitCheckpointSavedEvent>;

using MCMCObserver = std::function<void(const MCMCEvent&)>;

}  // namespace gelex

#endif  // GELEX_INFRA_LOGGING_FIT_EVENT_H_
