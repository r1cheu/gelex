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

#ifndef GELEX_ENGINE_MCMC_H_
#define GELEX_ENGINE_MCMC_H_

#include <filesystem>
#include <optional>
#include <string>
#include <vector>

#include "gelex/algo/infer/params.h"
#include "gelex/infra/logging/fit_event.h"
#include "gelex/model/bayes/method.h"

namespace gelex
{

class PhenoPipe;
class GenoPipe;
class BayesModel;

namespace mcmc
{

class FitEngine
{
   public:
    struct Config
    {
        std::string bfile_prefix;
        bayes::BayesConfig method;

        int seed;
        mcmc::Params mcmc_params;

        std::optional<std::vector<double>> pi;
        std::optional<std::vector<double>> dpi;
        std::optional<std::vector<double>> multiplier;
        std::optional<std::vector<double>> dmultiplier;

        double positive_prob{0.5};
        std::string out_prefix;
        std::optional<std::filesystem::path> resume_path;
    };

    explicit FitEngine(Config config);
    auto run(
        PhenoPipe&& pheno,
        GenoPipe&& geno,
        const FitObserver& observer = {}) -> void;

   private:
    Config config_;
};

// NOLINTBEGIN(cppcoreguidelines-avoid-const-or-ref-data-members)
class ConfigValidator
{
   public:
    explicit ConfigValidator(const FitEngine::Config& config) : config_{config}
    {
    }

    auto validate() const -> void;

   private:
    auto check_method() const -> void;
    auto check_mcmc_params() const -> void;
    auto check_mixture_priors() const -> void;
    auto check_positive_prob() const -> void;

    const FitEngine::Config& config_;
};
// NOLINTEND(cppcoreguidelines-avoid-const-or-ref-data-members)

}  // namespace mcmc

}  // namespace gelex

#endif  // GELEX_ENGINE_MCMC_H_
