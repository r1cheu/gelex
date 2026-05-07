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
class BayesModel;

namespace mcmc
{
class Result;
}

struct GeneticSummary;

namespace bayes
{
struct GeneticEffect;
}  // namespace bayes

enum class GeneticMode : uint8_t;
}  // namespace gelex

namespace gelex::cli
{

class McmcReporter : public FitReporter
{
   public:
    McmcReporter();

    using FitReporter::on_event;

    auto on_event(const MCMCBannerEvent& event) const -> void;
    auto on_event(const MCMCConfigEvent& event) const -> void;
    auto on_event(const MCMCProgressEvent& event) -> void;
    auto on_event(const MCMCCompleteEvent& event) const -> void;
    auto on_event(const FitCheckpointSavedEvent& event) -> void;

    auto as_observer() -> MCMCObserver
    {
        return [this](const MCMCEvent& e)
        { std::visit([this](const auto& ev) { this->on_event(ev); }, e); };
    }

   private:
    auto print_fixed_summary(
        const mcmc::Result& result,
        std::ptrdiff_t samples_collected) const -> void;
    auto print_genetic_summary(
        const GeneticSummary* summary,
        const bayes::GeneticEffect* effect,
        GeneticMode type) const -> void;
    auto print_residual_summary(const mcmc::Result& result) const -> void;

    size_t iter_{0};
    ProgressBar bar_;
    bool init_progress_ = false;
    std::string stats_;
};

}  // namespace gelex::cli

#endif  // GELEX_CLI_MCMC_REPORTER_H_
