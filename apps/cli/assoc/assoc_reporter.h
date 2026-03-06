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

#ifndef GELEX_CLI_ASSOC_REPORTER_H_
#define GELEX_CLI_ASSOC_REPORTER_H_

#include <cstddef>
#include <memory>

#include "gelex/infra/detail/indicator.h"
#include "gelex/infra/logging/assoc_event.h"
#include "gelex/infra/utils/utils.h"

namespace gelex
{
struct AssocConfigLoadedEvent;
struct AssocRemlStartedEvent;
struct AssocScanSummaryEvent;
struct AssocScanProgressEvent;
struct AssocLocoPhaseEvent;
struct AssocLocoRemlSummaryEvent;
struct AssocCompleteEvent;
}  // namespace gelex

namespace spdlog
{
class logger;
}

namespace gelex::cli
{

class AssocReporter
{
   public:
    AssocReporter();

    auto on_event(const AssocConfigLoadedEvent& event) const -> void;
    auto on_event(const AssocRemlStartedEvent& event) const -> void;
    auto on_event(const AssocScanSummaryEvent& event) -> void;
    auto on_event(const AssocScanProgressEvent& event) -> void;
    auto on_event(const AssocLocoPhaseEvent& event) -> void;
    auto on_event(const AssocLocoRemlSummaryEvent& event) -> void;
    auto on_event(const AssocCompleteEvent& event) -> void;

    auto as_observer() -> AssocObserver
    {
        return [this](const AssocEvent& e)
        { std::visit([this](const auto& ev) { this->on_event(ev); }, e); };
    }

   private:
    std::shared_ptr<spdlog::logger> logger_;
    size_t progress_{0};
    detail::ProgressBar bar_;
    bool bar_active_ = false;
    SmoothEtaCalculator eta_;
};

}  // namespace gelex::cli

#endif  // GELEX_CLI_ASSOC_REPORTER_H_
