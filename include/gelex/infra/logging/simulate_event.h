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
#ifndef GELEX_INFRA_LOGGING_SIMULATE_EVENT_H_
#define GELEX_INFRA_LOGGING_SIMULATE_EVENT_H_

#include <functional>
#include <optional>
#include <string>
#include <variant>

namespace gelex
{

struct HeritabilityGeneratedEvent
{
    double additive{};
    std::optional<double> dominance;
};

struct SimulateProgressEvent
{
    size_t total;
    size_t current;
    bool done;
};

struct SimulateConfigLoadedEvent
{
    double intercept{};
    double add_heritability{};
    std::optional<double> dom_heritability;
    int seed{};
};

struct OutputsWrittenEvent
{
    std::string phenotype_path;
    std::string snp_effect_path;
};

using SimulateEvent = std::variant<
    SimulateConfigLoadedEvent,
    SimulateProgressEvent,
    HeritabilityGeneratedEvent,
    OutputsWrittenEvent>;

using SimulateObserver = std::function<void(const SimulateEvent& event)>;

}  // namespace gelex

#endif  // GELEX_INFRA_LOGGING_SIMULATE_EVENT_H_
