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

#include "gelex/engine/mcmc.h"

#include <optional>
#include <string>
#include <utility>

#include <fmt/format.h>

#include "gelex/algo/mcmc/solver.h"
#include "gelex/bayes/prior.h"
#include "gelex/exception.h"
#include "gelex/infra/logging/fit_event.h"
#include "gelex/infra/logging/notify.h"
#include "gelex/io/mcmc.h"

namespace gelex
{

namespace mcmc
{

Engine::Engine(Config config) : config_(std::move(config)) {}

auto ConfigValidator::validate() const -> void
{
    check_mcmc_params();
}

auto ConfigValidator::check_mcmc_params() const -> void
{
    const auto& p = config_.mcmc_params;
    if (p.n_iters <= 0)
    {
        throw GelexException(
            fmt::format("n_iters must be positive, got {}", p.n_iters));
    }
    if (p.n_burn_in < 0 || p.n_burn_in >= p.n_iters)
    {
        throw GelexException(
            fmt::format(
                "n_burn_in must satisfy 0 <= n_burn_in < n_iters, got {} "
                "(n_iters={})",
                p.n_burn_in,
                p.n_iters));
    }
    if (p.n_thin <= 0)
    {
        throw GelexException(
            fmt::format("n_thin must be positive, got {}", p.n_thin));
    }
    if ((p.n_iters - p.n_burn_in) % p.n_thin != 0)
    {
        throw GelexException(
            fmt::format(
                "n_thin ({}) must divide n_iters - n_burn_in ({})",
                p.n_thin,
                p.n_iters - p.n_burn_in));
    }
    if (p.checkpoint_step < 0)
    {
        throw GelexException(
            fmt::format(
                "checkpoint_step must be non-negative, got {}",
                p.checkpoint_step));
    }
}

auto Engine::run(
    const BayesModel& model,
    bayes::BayesPrior prior,
    const MCMCObserver& observer) -> void
{
    notify(observer, FitPriorSetEvent{&prior});

    std::optional<std::string> checkpoint_prefix;
    if (config_.mcmc_params.checkpoint_step > 0)
    {
        checkpoint_prefix = config_.out_prefix;
    }
    auto solver = mcmc::Solver{
        config_.mcmc_params,
        config_.out_prefix + ".draws",
        std::move(checkpoint_prefix)};
    auto result = [&]() -> mcmc::Result
    {
        if (config_.from_checkpoint_path)
        {
            return solver.run_from(
                model,
                std::move(prior),
                *config_.from_checkpoint_path,
                observer);
        }
        return solver.run(model, prior, config_.seed, observer);
    }();
    mcmc::write_params(result, config_.out_prefix);
    mcmc::write_summary(result, config_.out_prefix);
    mcmc::write_snp_eff(
        result, model, config_.bfile_prefix + ".bim", config_.out_prefix);
    notify(observer, FitResultsSavedEvent{config_.out_prefix});
}

}  // namespace mcmc

}  // namespace gelex
