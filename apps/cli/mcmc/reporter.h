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
#include <string_view>

#include "gelex/bayes/genotype/progress.h"
#include "gelex/bayes/mcmc_progress.h"
#include "gelex/genetic_mode.h"

#include "cli/progress_bar.h"
#include "cli/timer.h"

namespace gelex
{
namespace bayes
{
class GeneticDesign;
}  // namespace bayes

class BayesModel;
}  // namespace gelex

namespace cli
{

class GenoReporter
{
   public:
    GenoReporter();

    auto show_total(int64_t num_variants) const -> void;
    auto show_loaded(const gelex::bayes::GeneticDesign& design) const -> void;
    auto on_event(const gelex::GenotypeProgressEvent& event) -> void;

   private:
    auto show_loaded_mode(
        gelex::GeneticMode mode,
        int64_t num_snps,
        int64_t invalid_snps) const -> void;

    std::size_t progress_{0};
    std::size_t total_{0};
    cli::ProgressBar bar_;
    bool bar_active_ = false;
    gelex::SmoothEtaCalculator eta_;
};

class McmcReporter
{
   public:
    McmcReporter() = default;

    auto show_dataset_summary(
        const gelex::BayesModel& model,
        std::string_view pheno_name) -> void;
    auto on_event(const gelex::MCMCProgressEvent& event) -> void;

   private:
    std::size_t iter_{0};
    cli::ProgressBar bar_;
    bool init_progress_ = false;
};

}  // namespace cli

#endif  // APPS_CLI_MCMC_REPORTER_H_
