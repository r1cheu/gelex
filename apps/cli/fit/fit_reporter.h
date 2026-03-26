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

#ifndef GELEX_CLI_FIT_REPORTER_H_
#define GELEX_CLI_FIT_REPORTER_H_

#include <cstddef>
#include <memory>
#include <string>
#include <string_view>

#include <Eigen/Core>

#include "gelex/infra/logging/fit_event.h"
#include "gelex/infra/logging/progress_bar.h"

namespace gelex
{
class BayesModel;
class MCMCResult;
struct GeneticSummary;
struct PosteriorSummary;

namespace detail
{
struct ScaledInvChiSqParams;
}  // namespace detail

namespace bayes
{
struct GeneticEffect;
struct GeneticPrior;
struct VariancePrior;
}  // namespace bayes

enum class GeneticKind : uint8_t;
}  // namespace gelex

namespace spdlog
{
class logger;
}

namespace gelex::cli
{

class FitReporter
{
   public:
    FitReporter();

    auto on_event(const FitConfigLoadedEvent& event) const -> void;
    auto on_event(const FitPriorSetEvent& event) const -> void;
    auto on_event(const FitMcmcProgressEvent& event) -> void;
    auto on_event(const FitMcmcCompleteEvent& event) const -> void;
    auto on_event(const FitResultsSavedEvent& event) const -> void;
    auto on_event(const FitCheckpointSavedEvent& event) -> void;

    auto as_observer() -> FitObserver
    {
        return [this](const FitEvent& e)
        { std::visit([this](const auto& ev) { this->on_event(ev); }, e); };
    }

   private:
    auto print_random_prior(const bayes::VariancePrior& prior) const -> void;
    auto print_genetic_prior(const bayes::GeneticPrior& prior) const -> void;
    auto print_residual_prior(const bayes::VariancePrior& prior) const -> void;

    auto print_fixed_summary(
        const MCMCResult& result,
        std::ptrdiff_t samples_collected) const -> void;
    auto print_genetic_summary(
        const GeneticSummary* summary,
        const bayes::GeneticEffect* effect,
        GeneticKind type) const -> void;
    auto print_residual_summary(const MCMCResult& result) const -> void;

    auto print_variance_prior(
        const detail::ScaledInvChiSqParams& prior,
        double init_variance) const -> void;
    auto print_summary_row(
        std::string_view name,
        const PosteriorSummary& summary,
        Eigen::Index index = 0) const -> void;

    std::shared_ptr<spdlog::logger> logger_;
    size_t iter_{0};
    ProgressBar bar_;
    bool init_progress_ = false;
    std::string stats_;
};

}  // namespace gelex::cli

#endif  // GELEX_CLI_FIT_REPORTER_H_
