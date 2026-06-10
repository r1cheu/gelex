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

#ifndef APPS_CLI_GENO_REPORTER_H_
#define APPS_CLI_GENO_REPORTER_H_

#include <variant>

#include "gelex/infra/logging/geno_event.h"
#include "gelex/infra/logging/progress_bar.h"

namespace cli
{

class GenoReporter
{
   public:
    GenoReporter();

    auto on_event(const gelex::GenotypeLoadedEvent& event) const -> void;
    auto on_event(const gelex::GenotypeProgressEvent& event) -> void;

    auto as_observer() -> gelex::GenoObserver
    {
        return [this](const gelex::GenoEvent& e)
        { std::visit([this](const auto& ev) { this->on_event(ev); }, e); };
    }

   private:
    gelex::ProgressInfo progress_info_;
    bool init_progress_ = false;
};

}  // namespace cli

#endif  // APPS_CLI_GENO_REPORTER_H_
