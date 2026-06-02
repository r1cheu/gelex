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

#include "gelex/algo/infer/mcmc/solver.h"

#include <cstddef>
#include <random>
#include <utility>

#include "gelex/algo/infer/mcmc/chain.h"
#include "gelex/algo/infer/mcmc/records.h"
#include "gelex/bayes/state.h"
#include "gelex/exception.h"
#include "gelex/infra/logging/notify.h"

namespace gelex::mcmc
{

Solver::Solver(
    mcmc::Params params,
    std::string sample_prefix,
    std::optional<std::string> checkpoint_prefix)
    : params_(params),
      sample_prefix_(std::move(sample_prefix)),
      checkpoint_prefix_(std::move(checkpoint_prefix))
{
}

auto Solver::run(
    const BayesModel& model,
    bayes::BayesPrior prior,
    Eigen::Index seed,
    const MCMCObserver& observer) -> mcmc::Result
{
    if (checkpoint_prefix_)
    {
        throw GelexException(
            "MCMC checkpoint writer is not implemented after Bayes prior/state "
            "cleanup");
    }
    if (!sample_prefix_.empty())
    {
        throw GelexException(
            "MCMC sample output is legacy and will be removed");
    }

    auto state = BayesState{model, prior};
    auto rng = std::mt19937_64{static_cast<std::mt19937_64::result_type>(seed)};
    auto records = mcmc::Records{};
    auto chain = Chain::make(model, prior, state, rng);

    for (Eigen::Index iter = 0; iter < params_.n_iters; ++iter)
    {
        chain.step();

        notify(
            observer,
            MCMCProgressEvent{
                .current = static_cast<size_t>(iter + 1),
                .total = static_cast<size_t>(params_.n_iters),
                .state = &state,
            });

        if (iter >= params_.n_burn_in
            && (iter + 1 - params_.n_burn_in) % params_.n_thin == 0)
        {
            records.store(model, prior, state);
        }
    }

    notify(
        observer,
        MCMCProgressEvent{
            .current = static_cast<size_t>(params_.n_iters),
            .total = static_cast<size_t>(params_.n_iters),
            .done = true,
        });

    auto result = mcmc::Result{std::move(records), model, params_.n_records()};
    notify(observer, MCMCCompleteEvent{&result, &model, params_.n_records()});
    return result;
}

auto Solver::resume(
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

}  // namespace gelex::mcmc
