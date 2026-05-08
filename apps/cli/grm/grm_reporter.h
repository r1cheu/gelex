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

#ifndef GELEX_CLI_GRM_REPORTER_H_
#define GELEX_CLI_GRM_REPORTER_H_

#include <cstddef>

#include "gelex/infra/logging/grm_event.h"
#include "gelex/infra/logging/progress_bar.h"
#include "gelex/infra/logging/timer.h"

namespace gelex::cli
{

class GrmReporter
{
   public:
    GrmReporter();

    auto on_event(const GrmBannerEvent& event) const -> void;
    auto on_event(const GrmConfigLoadedEvent& event) const -> void;
    auto on_event(const GrmDataLoadedEvent& event) const -> void;
    auto on_event(const GrmComputeStartedEvent& event) -> void;
    auto on_event(const GrmProgressEvent& event) -> void;
    auto on_event(const GrmFilesWrittenEvent& event) const -> void;

    auto as_observer() -> GrmObserver
    {
        return [this](const GrmEvent& e)
        { std::visit([this](const auto& ev) { this->on_event(ev); }, e); };
    }

   private:
    size_t progress_{0};
    ProgressBar bar_;
    bool bar_active_ = false;
    SmoothEtaCalculator eta_;
    size_t global_total_ = 0;
    size_t accumulated_base_ = 0;
};

}  // namespace gelex::cli

#endif  // GELEX_CLI_GRM_REPORTER_H_
