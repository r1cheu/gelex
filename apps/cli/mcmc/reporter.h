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

#include "cli/progress.h"

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
    explicit GenoReporter(std::size_t total);

    static auto show_total(int64_t num_variants) -> void;
    auto show_loaded(const gelex::bayes::GeneticDesign& design) const -> void;
    auto on_event(const gelex::GenotypeProgressEvent& event) -> void;

   private:
    auto show_loaded_mode(
        gelex::GeneticMode mode,
        int64_t num_snps,
        int64_t invalid_snps) const -> void;

    cli::Progress progress_;
    decltype(cli::make_rate()) estimate_rate_;
    decltype(cli::make_eta(std::size_t{})) estimate_eta_;
};

class McmcReporter
{
   public:
    McmcReporter(std::size_t total, std::size_t burn_in);

    static auto show_dataset_summary(
        const gelex::BayesModel& model,
        std::string_view pheno_name) -> void;
    static auto show_sampling_header() -> void;
    auto on_event(const gelex::MCMCProgressEvent& event) -> void;

   private:
    std::size_t burn_in_;
    cli::Progress progress_;
    decltype(cli::make_rate()) estimate_rate_;
    decltype(cli::make_eta(std::size_t{})) estimate_eta_;
};

}  // namespace cli

#endif  // APPS_CLI_MCMC_REPORTER_H_
