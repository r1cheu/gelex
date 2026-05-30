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

#ifndef GELEX_ALGO_INFER_MCMC_SOLVER_H_
#define GELEX_ALGO_INFER_MCMC_SOLVER_H_

#include <filesystem>
#include <optional>
#include <string>
#include <utility>

#include <Eigen/Core>

#include "gelex/algo/infer/mcmc/result.h"
#include "gelex/algo/infer/params.h"
#include "gelex/bayes/model.h"
#include "gelex/bayes/prior.h"
#include "gelex/exception.h"
#include "gelex/infra/logging/fit_event.h"

namespace gelex
{

namespace mcmc
{

template <typename ChainFactory>
class Solver
{
   public:
    explicit Solver(
        mcmc::Params params,
        ChainFactory make_chain,
        std::string sample_prefix = {},
        std::optional<std::string> checkpoint_prefix = std::nullopt);

    auto run(
        const BayesModel& model,
        bayes::BayesPrior prior,
        Eigen::Index seed = 42,
        const MCMCObserver& observer = {}) -> mcmc::Result;

    auto resume(
        const BayesModel& model,
        bayes::BayesPrior prior,
        const std::filesystem::path& checkpoint_path,
        const MCMCObserver& observer = {}) -> mcmc::Result;

   private:
    [[no_unique_address]] ChainFactory make_chain_;
    mcmc::Params params_;
    std::string sample_prefix_;
    std::optional<std::string> checkpoint_prefix_;
};

template <typename ChainFactory>
Solver<ChainFactory>::Solver(
    mcmc::Params params,
    ChainFactory make_chain,
    std::string sample_prefix,
    std::optional<std::string> checkpoint_prefix)
    : make_chain_(std::move(make_chain)),
      params_(params),
      sample_prefix_(std::move(sample_prefix)),
      checkpoint_prefix_(std::move(checkpoint_prefix))
{
}

template <typename ChainFactory>
auto Solver<ChainFactory>::run(
    const BayesModel& model,
    bayes::BayesPrior prior,
    Eigen::Index seed,
    const MCMCObserver& observer) -> mcmc::Result
{
    static_cast<void>(model);
    static_cast<void>(prior);
    static_cast<void>(seed);
    static_cast<void>(observer);
    throw GelexException(
        "MCMC solver is not implemented after Bayes prior/state cleanup");
}

template <typename ChainFactory>
auto Solver<ChainFactory>::resume(
    const BayesModel& model,
    bayes::BayesPrior prior,
    const std::filesystem::path& checkpoint_path,
    const MCMCObserver& observer) -> mcmc::Result
{
    static_cast<void>(model);
    static_cast<void>(prior);
    static_cast<void>(checkpoint_path);
    static_cast<void>(observer);
    throw GelexException(
        "MCMC solver resume is not implemented after Bayes prior/state "
        "cleanup");
}

}  // namespace mcmc

}  // namespace gelex

#endif  // GELEX_ALGO_INFER_MCMC_SOLVER_H_
