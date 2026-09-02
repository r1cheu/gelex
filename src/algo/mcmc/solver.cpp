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

#include "gelex/algo/mcmc/solver.h"

#include <cstddef>
#include <fmt/format.h>
#include <random>
#include <utility>

#include "gelex/algo/mcmc/chain.h"
#include "gelex/algo/mcmc/records.h"
#include "gelex/bayes/state.h"
#include "gelex/exception.h"
#include "gelex/infra/logging/notify.h"
#include "gelex/io/mcmc_checkpoint.h"

namespace gelex
{

Solver::Solver(
    MCMCParams params,
    std::string draws_path,
    std::optional<std::string> checkpoint_prefix)
    : params_(params),
      draws_path_(std::move(draws_path)),
      checkpoint_prefix_(std::move(checkpoint_prefix))
{
    if (params_.n_iters <= 0)
    {
        throw GelexException(
            fmt::format("n_iters must be positive, got {}", params_.n_iters));
    }
    if (params_.n_burn_in < 0 || params_.n_burn_in >= params_.n_iters)
    {
        throw GelexException(
            fmt::format(
                "n_burn_in must satisfy 0 <= n_burn_in < n_iters, got {} "
                "(n_iters={})",
                params_.n_burn_in,
                params_.n_iters));
    }
    if (params_.n_thin <= 0)
    {
        throw GelexException(
            fmt::format("n_thin must be positive, got {}", params_.n_thin));
    }
    if ((params_.n_iters - params_.n_burn_in) % params_.n_thin != 0)
    {
        throw GelexException(
            fmt::format(
                "n_thin ({}) must divide n_iters - n_burn_in ({})",
                params_.n_thin,
                params_.n_iters - params_.n_burn_in));
    }
    if (params_.checkpoint_step < 0)
    {
        throw GelexException(
            fmt::format(
                "checkpoint_step must be non-negative, got {}",
                params_.checkpoint_step));
    }
}

auto Solver::run(
    const BayesModel& model,
    const bayes::BayesPrior& prior,
    Eigen::Index seed,
    const MCMCObserver& observer) -> Result
{
    auto state = BayesState{model, prior};
    auto rng = std::mt19937_64{static_cast<std::mt19937_64::result_type>(seed)};
    return run_iterations(model, prior, state, rng, observer);
}

auto Solver::run_from(
    const BayesModel& model,
    const bayes::BayesPrior& prior,
    const std::filesystem::path& checkpoint_path,
    const MCMCObserver& observer) -> Result
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
    const MCMCObserver& observer) -> Result
{
    auto records = Records{params_.n_records(), draws_path_};
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
            notify(observer, MCMCCheckpointSavedEvent{});
        }
    }

    notify(
        observer,
        MCMCProgressEvent{
            .current = static_cast<size_t>(params_.n_iters),
            .total = static_cast<size_t>(params_.n_iters),
            .done = true,
        });

    return Result{std::move(records), model, params_.n_records()};
}

}  // namespace gelex
