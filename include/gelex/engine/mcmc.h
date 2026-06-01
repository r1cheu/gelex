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
#include "gelex/bayes/prior.h"
#include "gelex/bayes/recipe.h"
#include "gelex/infra/logging/fit_event.h"

namespace gelex
{

class BayesModel;

namespace mcmc
{

class Engine
{
   public:
    struct Config
    {
        std::string bfile_prefix;
        bayes::BayesRecipeScheme recipe_scheme;
        bayes::BayesRecipeConfig recipe_config;

        int seed;
        mcmc::Params mcmc_params;

        std::string out_prefix;
        std::optional<std::filesystem::path> resume_path;
    };

    explicit Engine(Config config);
    auto run(
        const BayesModel& model,
        bayes::BayesPrior prior,
        const MCMCObserver& observer = {}) -> void;

   private:
    Config config_;
};

// NOLINTBEGIN(cppcoreguidelines-avoid-const-or-ref-data-members)
class ConfigValidator
{
   public:
    explicit ConfigValidator(const Engine::Config& config) : config_{config} {}

    auto validate() const -> void;

   private:
    auto check_method() const -> void;
    auto check_mcmc_params() const -> void;

    const Engine::Config& config_;
};
// NOLINTEND(cppcoreguidelines-avoid-const-or-ref-data-members)

}  // namespace mcmc

}  // namespace gelex

#endif  // GELEX_ENGINE_MCMC_H_
