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
#include <optional>
#include <string>
#include <variant>

#include "gelex/types/effects.h"
#include "gelex/types/genetic_effect_type.h"

namespace gelex
{

class BayesModel;
class MCMCResult;

struct FitConfigLoadedEvent
{
    gelex::BayesAlphabet method;
    ModelType model_type;
    int n_iters;
    int n_burnin;
    int seed;
};

struct FitModelReadyEvent
{
    const BayesModel* model;
};

struct FitMcmcProgressEvent
{
    size_t current{};
    size_t total{};
    bool done{};
    std::optional<double> h2;
    std::optional<double> d2;
    std::optional<double> sigma2_e;
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

using FitEvent = std::variant<
    FitConfigLoadedEvent,
    FitModelReadyEvent,
    FitMcmcProgressEvent,
    FitMcmcCompleteEvent,
    FitResultsSavedEvent>;

using FitObserver = std::function<void(const FitEvent&)>;

}  // namespace gelex

#endif  // GELEX_INFRA_LOGGING_FIT_EVENT_H_
