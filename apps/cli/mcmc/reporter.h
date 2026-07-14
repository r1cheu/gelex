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

#ifndef APPS_CLI_MCMC_REPORTER_H_
#define APPS_CLI_MCMC_REPORTER_H_

#include <cstddef>
#include <cstdint>
#include <string>
#include <string_view>

#include "gelex/bayes/prior.h"
#include "gelex/infra/logging/fit_event.h"
#include "gelex/infra/logging/geno_event.h"
#include "gelex/types/genetic_mode.h"

#include "cli/progress_bar.h"

namespace gelex
{
namespace bayes
{
class RandomPrior;
class ResidualPrior;
class ScaledInvChiSqPrior;
}  // namespace bayes

class BayesModel;
}  // namespace gelex

namespace cli
{

class GenoReporter
{
   public:
    GenoReporter();

    auto show_loaded(
        gelex::GeneticMode mode,
        int64_t num_snps,
        int64_t invalid_snps) const -> void;
    auto on_event(const gelex::GenotypeProgressEvent& event) -> void;

    auto as_observer() -> gelex::GenoObserver
    {
        return [this](const gelex::GenotypeProgressEvent& e)
        { this->on_event(e); };
    }

   private:
    cli::ProgressInfo progress_info_;
    bool init_progress_ = false;
};

class McmcReporter
{
   public:
    McmcReporter() = default;

    auto show_dataset_summary(
        const gelex::BayesModel& model,
        std::string_view pheno_name) -> void;
    auto show_prior(const gelex::bayes::BayesPrior& prior) -> void;
    auto show_complete(std::ptrdiff_t samples_collected) -> void;
    auto on_event(const gelex::MCMCProgressEvent& event) -> void;
    auto on_event(const gelex::MCMCCheckpointSavedEvent& event) -> void;

    auto as_observer() -> gelex::MCMCObserver
    {
        return [this](const gelex::MCMCEvent& e)
        { std::visit([this](const auto& ev) { this->on_event(ev); }, e); };
    }

   private:
    static auto print_random_prior(const gelex::bayes::RandomPrior& prior)
        -> void;
    static auto print_genetic_prior(const gelex::bayes::GeneticPrior& prior)
        -> void;
    static auto print_residual_prior(const gelex::bayes::ResidualPrior& prior)
        -> void;
    static auto print_variance_prior(
        const gelex::bayes::ScaledInvChiSqPrior& prior,
        double init_variance) -> void;

    size_t iter_{0};
    cli::ProgressBar bar_;
    bool init_progress_ = false;
    std::string stats_;
    const gelex::bayes::BayesPrior* prior_ = nullptr;
};

}  // namespace cli

#endif  // APPS_CLI_MCMC_REPORTER_H_
