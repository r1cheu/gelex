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

#include <functional>
#include <optional>
#include <string>
#include <variant>

#include "gelex/types/effects.h"

namespace gelex
{

struct FitConfigLoadedEvent
{
    gelex::BayesAlphabet method;
    bool use_dominance;
    int n_iters;
    int n_burnin;
    int seed;
    int threads;
};

struct FitMcmcProgressEvent
{
    size_t current{};
    size_t total{};
    bool done{};
    std::optional<double> h2;
    std::optional<double> h2_dom;
    std::optional<double> sigma2_e;
};

struct FitResultsSavedEvent
{
    std::string out_prefix;
};

using FitEvent = std::
    variant<FitConfigLoadedEvent, FitMcmcProgressEvent, FitResultsSavedEvent>;

using FitObserver = std::function<void(const FitEvent&)>;

}  // namespace gelex

#endif  // GELEX_INFRA_LOGGING_FIT_EVENT_H_
