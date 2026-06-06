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

#ifndef GELEX_CLI_MCMC_REPORTER_H_
#define GELEX_CLI_MCMC_REPORTER_H_

#include <cstddef>
#include <string>

#include "cli/fit_reporter.h"
#include "gelex/infra/logging/fit_event.h"
#include "gelex/infra/logging/progress_bar.h"

namespace gelex
{
namespace cli
{

class McmcReporter : public FitReporter
{
   public:
    McmcReporter() = default;

    using FitReporter::on_event;

    auto on_event(const MCMCProgressEvent& event) -> void;
    auto on_event(const FitCheckpointSavedEvent& event) -> void;

    auto as_observer() -> MCMCObserver
    {
        return [this](const MCMCEvent& e)
        { std::visit([this](const auto& ev) { this->on_event(ev); }, e); };
    }

   private:
    auto on_event(const MCMCBannerEvent& event) -> void;
    auto on_event(const MCMCConfigEvent& event) -> void;
    auto on_event(const FitPriorSetEvent& event) -> void;
    auto on_event(const MCMCCompleteEvent& event) -> void;
    auto on_event(const FitResultsSavedEvent& event) -> void;

    size_t iter_{0};
    ProgressBar bar_;
    bool init_progress_ = false;
    std::string stats_;
    const bayes::BayesPrior* prior_ = nullptr;
};

}  // namespace cli

}  // namespace gelex

#endif  // GELEX_CLI_MCMC_REPORTER_H_
