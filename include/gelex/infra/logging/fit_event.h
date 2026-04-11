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

#include "gelex/types/bayes_method.h"
#include "gelex/types/genetic_effect_type.h"

namespace gelex
{

class BayesModel;
class BayesState;
class MCMCResult;

namespace bayes
{
class Priors;
}

struct FitConfigLoadedEvent
{
    gelex::BayesMethodConfig method;
    GeneticMode model_type{};
    int n_iters{};
    int n_burn_in{};
    int seed{};
};

struct FitPriorSetEvent
{
    const bayes::Priors* priors{};
};

struct FitMcmcProgressEvent
{
    size_t current{};
    size_t total{};
    bool done{};
    const BayesState* state{};
};

struct FitMcmcCompleteEvent
{
    const MCMCResult* result;
    const BayesModel* model;
    std::ptrdiff_t samples_collected;
};

struct FitResultsSavedEvent
{
    std::string out_prefix;
};

struct FitCheckpointSavedEvent
{
};

using FitEvent = std::variant<
    FitConfigLoadedEvent,
    FitPriorSetEvent,
    FitMcmcProgressEvent,
    FitMcmcCompleteEvent,
    FitResultsSavedEvent,
    FitCheckpointSavedEvent>;

using FitObserver = std::function<void(const FitEvent&)>;

}  // namespace gelex

#endif  // GELEX_INFRA_LOGGING_FIT_EVENT_H_
