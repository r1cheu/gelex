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

#ifndef APPS_CLI_MCMC_PROGRESS_H_
#define APPS_CLI_MCMC_PROGRESS_H_

#include <cstddef>

#include "cli/progress.h"

namespace cli
{

class GenotypeProgress
{
   public:
    explicit GenotypeProgress(std::size_t total);
    auto operator()(std::size_t current) -> void;
    auto finish() -> void;

   private:
    cli::Progress progress_;
    decltype(cli::make_rate()) estimate_rate_;
    decltype(cli::make_eta(std::size_t{})) estimate_eta_;
};

class McmcProgress
{
   public:
    McmcProgress(std::size_t total, int burn_in);
    auto operator()(std::size_t current) -> void;
    auto finish() -> void;

   private:
    std::size_t burn_in_;
    cli::Progress progress_;
    decltype(cli::make_rate()) estimate_rate_;
    decltype(cli::make_eta(std::size_t{})) estimate_eta_;
};

}  // namespace cli

#endif  // APPS_CLI_MCMC_PROGRESS_H_
