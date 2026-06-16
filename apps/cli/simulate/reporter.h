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

#ifndef APPS_CLI_SIMULATE_REPORTER_H_
#define APPS_CLI_SIMULATE_REPORTER_H_

#include <optional>

#include "gelex/infra/logging/progress_bar.h"
#include "gelex/infra/logging/simulate_event.h"

namespace cli
{

class SimulatorReporter
{
   public:
    SimulatorReporter();

    auto show_banner() const -> void;
    auto show_config(
        std::optional<double> add_heritability,
        std::optional<double> dom_heritability,
        int seed) const -> void;
    auto on_event(const gelex::SimulateProgressEvent& event) -> void;
    auto show_variance_summary(
        std::optional<double> realized_h2,
        std::optional<double> realized_d2) const -> void;

    auto as_observer() -> gelex::SimulateObserver
    {
        return [this](const gelex::SimulateProgressEvent& e)
        { this->on_event(e); };
    }

   private:
    gelex::ProgressInfo info_;
    bool init_progress_ = false;
};

}  // namespace cli

#endif  // APPS_CLI_SIMULATE_REPORTER_H_
