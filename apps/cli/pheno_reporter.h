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

#ifndef APPS_CLI_PHENO_REPORTER_H_
#define APPS_CLI_PHENO_REPORTER_H_

#include <variant>

#include "gelex/infra/logging/pheno_event.h"

namespace cli
{

class PhenoReporter
{
   public:
    auto on_event(const gelex::PhenotypeLoadedEvent& event) const -> void;
    auto on_event(const gelex::CovariatesLoadedEvent& event) const -> void;

    auto as_observer() -> gelex::PhenoObserver
    {
        return [this](const gelex::PhenoEvent& e)
        { std::visit([this](const auto& ev) { this->on_event(ev); }, e); };
    }
};

}  // namespace cli

#endif  // APPS_CLI_PHENO_REPORTER_H_
