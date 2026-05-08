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

#ifndef GELEX_CLI_DATASET_REPORTER_H_
#define GELEX_CLI_DATASET_REPORTER_H_

#include <variant>

#include "gelex/infra/logging/dataset_event.h"

namespace gelex::cli
{

class DatasetReporter
{
   public:
    auto on_event(const DatasetSectionEvent& event) const -> void;
    auto on_event(const IntersectionEvent& event) const -> void;

    auto as_observer() -> DatasetObserver
    {
        return [this](const DatasetEvent& e)
        { std::visit([this](const auto& ev) { this->on_event(ev); }, e); };
    }
};

}  // namespace gelex::cli

#endif  // GELEX_CLI_DATASET_REPORTER_H_
