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

#ifndef GELEX_CLI_VI_REPORTER_H_
#define GELEX_CLI_VI_REPORTER_H_

#include "cli/fit_reporter.h"
#include "gelex/infra/logging/fit_event.h"
#include "gelex/infra/logging/progress_bar.h"

namespace gelex::cli
{

class ViReporter : public FitReporter
{
   public:
    ViReporter();

    using FitReporter::on_event;

    auto on_event(const VIBannerEvent& event) const -> void;
    auto on_event(const VIConfigEvent& event) const -> void;
    auto on_event(const VIProgressEvent& event) -> void;
    auto on_event(const VICompleteEvent& event) const -> void;

    auto as_observer() -> VIObserver
    {
        return [this](const VIEvent& e)
        { std::visit([this](const auto& ev) { this->on_event(ev); }, e); };
    }

   private:
    ProgressInfo cavi_info_;
    bool init_progress_ = false;
};

}  // namespace gelex::cli

#endif  // GELEX_CLI_VI_REPORTER_H_
