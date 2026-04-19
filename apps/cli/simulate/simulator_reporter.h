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

#ifndef GELEX_CLI_SIMULATOR_REPORTER_H_
#define GELEX_CLI_SIMULATOR_REPORTER_H_

#include <memory>

#include "gelex/infra/logging/progress_bar.h"
#include "gelex/infra/logging/simulate_event.h"

namespace spdlog
{
class logger;
}
namespace gelex::cli
{

class SimulatorReporter
{
   public:
    SimulatorReporter();

    auto on_event(const SimulateBannerEvent& event) const -> void;
    auto on_event(const SimulateConfigLoadedEvent& event) const -> void;
    auto on_event(const SimulateProgressEvent& event) -> void;

    auto as_observer() -> SimulateObserver
    {
        return [this](const SimulateEvent& e)
        { std::visit([this](const auto& ev) { this->on_event(ev); }, e); };
    }

   private:
    std::shared_ptr<spdlog::logger> logger_;
    ProgressInfo info_;
    bool init_progress_ = false;
};

}  // namespace gelex::cli

#endif  // GELEX_CLI_SIMULATOR_REPORTER_H_
