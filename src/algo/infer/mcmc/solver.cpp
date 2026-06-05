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
#include "gelex/infra/logging/notify.h"
#include "gelex/io/mcmc/checkpoint_reader.h"
#include "gelex/io/mcmc/checkpoint_writer.h"

namespace gelex::mcmc
{

Solver::Solver(
    mcmc::Params params,
    std::string draws_path,
    std::optional<std::string> checkpoint_prefix)
    : params_(params),
      draws_path_(std::move(draws_path)),
      checkpoint_prefix_(std::move(checkpoint_prefix))
{
}

auto Solver::run(
    const BayesModel& model,
    const bayes::BayesPrior& prior,
    Eigen::Index seed,
    const MCMCObserver& observer) -> mcmc::Result
{
    auto state = BayesState{model, prior};
    auto rng = std::mt19937_64{static_cast<std::mt19937_64::result_type>(seed)};
    return run_iterations(model, prior, state, rng, observer);
}

auto Solver::run_from(
    const BayesModel& model,
    bayes::BayesPrior prior,
    const std::filesystem::path& checkpoint_path,
    const MCMCObserver& observer) -> mcmc::Result
{
    auto state = BayesState{model, prior};
    auto rng = read_checkpoint(checkpoint_path, state);
    return run_iterations(model, prior, state, rng, observer);
}

auto Solver::run_iterations(
    const BayesModel& model,
    const bayes::BayesPrior& prior,
    BayesState& state,
    std::mt19937_64& rng,
    const MCMCObserver& observer) -> mcmc::Result
{
    auto records = mcmc::Records{params_.n_records(), draws_path_};
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
            records.store(model, state);
        }

        if (checkpoint_prefix_ && params_.checkpoint_step > 0
            && (iter + 1) % params_.checkpoint_step == 0)
        {
            write_checkpoint(state, rng, *checkpoint_prefix_);
            notify(observer, FitCheckpointSavedEvent{});
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

}  // namespace gelex::mcmc
