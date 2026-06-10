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

#ifndef APPS_CLI_ASSOC_REPORTER_H_
#define APPS_CLI_ASSOC_REPORTER_H_

#include <cstddef>

#include "gelex/infra/logging/assoc_event.h"
#include "gelex/infra/logging/progress_bar.h"
#include "gelex/infra/logging/timer.h"

namespace cli
{

class AssocReporter
{
   public:
    AssocReporter();

    auto on_event(const gelex::AssocBannerEvent& event) const -> void;
    auto on_event(const gelex::AssocConfigLoadedEvent& event) const -> void;
    auto on_event(const gelex::AssocRemlStartedEvent& event) const -> void;
    auto on_event(const gelex::AssocScanSummaryEvent& event) -> void;
    auto on_event(const gelex::AssocScanProgressEvent& event) -> void;
    auto on_event(const gelex::AssocLocoPhaseEvent& event) -> void;
    auto on_event(const gelex::AssocLocoRemlSummaryEvent& event) -> void;
    auto on_event(const gelex::AssocCompleteEvent& event) -> void;

    auto as_observer() -> gelex::AssocObserver
    {
        return [this](const gelex::AssocEvent& e)
        { std::visit([this](const auto& ev) { this->on_event(ev); }, e); };
    }

   private:
    size_t progress_{0};
    gelex::ProgressBar bar_;
    bool bar_active_ = false;
    gelex::SmoothEtaCalculator eta_;
};

}  // namespace cli

#endif  // APPS_CLI_ASSOC_REPORTER_H_
