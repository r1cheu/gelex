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
#include <string_view>
#include <utility>

#include <fmt/format.h>
#include <omp.h>
#include <Eigen/Core>

#include "gelex/algo/detail/posterior_calculator.h"
#include "gelex/algo/infer/mcmc/context.h"
#include "gelex/algo/infer/mcmc/result.h"
#include "gelex/algo/infer/mcmc/samples.h"
#include "gelex/algo/infer/params.h"
#include "gelex/data/genotype/genotype.h"
#include "gelex/exception.h"
#include "gelex/infra/detail/eigen_thread_guard.h"
#include "gelex/infra/logging/fit_event.h"
#include "gelex/infra/logging/notify.h"
#include "gelex/io/mcmc/checkpoint_reader.h"
#include "gelex/io/mcmc/checkpoint_writer.h"
#include "gelex/model/bayes/model.h"
#include "gelex/model/bayes/prior.h"

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
    static void validate_checkpoint(
        const mcmc::State& state,
        const BayesModel& model);

    template <typename Chain>
    auto run_impl(
        Chain& chain,
        mcmc::Samples& samples,
        mcmc::State& state,
        std::mt19937_64& rng,
        const MCMCObserver& observer) -> void;

    auto finalize(
        mcmc::Samples samples,
        const BayesModel& model,
        const MCMCObserver& observer) -> mcmc::Result;

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
    mcmc::State state{model, prior};
    mcmc::Samples samples(
        model, prior, state, sample_prefix_, params_.n_records());
    notify(observer, FitPriorSetEvent{.prior = &prior});
    std::mt19937_64 rng(seed);

    const infra::detail::EigenThreadGuard guard;
    omp_set_num_threads(1);

    mcmc::Context ctx{
        .model = model, .prior = prior, .state = state, .rng = rng};
    auto chain = make_chain_(ctx);
    run_impl(chain, samples, state, rng, observer);
    return finalize(std::move(samples), model, observer);
}

template <typename ChainFactory>
void Solver<ChainFactory>::validate_checkpoint(
    const mcmc::State& state,
    const BayesModel& model)
{
    auto check
        = [](Eigen::Index actual, Eigen::Index expected, std::string_view label)
    {
        if (actual != expected)
        {
            throw GelexException(
                fmt::format(
                    "checkpoint dimension mismatch for {}: expected {}, got {}",
                    label,
                    expected,
                    actual));
        }
    };

    check(state.fixed().coeffs.size(), model.fixed().X.cols(), "fixed.coeffs");

    if (state.random().size() != model.random().size())
    {
        throw GelexException(
            fmt::format(
                "checkpoint random design count mismatch: expected {}, got {}",
                model.random().size(),
                state.random().size()));
    }
    for (size_t i = 0; i < state.random().size(); ++i)
    {
        check(
            state.random()[i].coeffs.size(),
            model.random()[i].X.cols(),
            fmt::format("random[{}].coeffs", i));
    }

    for (const auto& gs : state.genetics())
    {
        const auto* design = model.genetic(gs.type);
        if (design == nullptr)
        {
            throw GelexException(
                fmt::format(
                    "checkpoint has genetic type {} not in model",
                    EffectType::from_genetic(gs.type)));
        }
        const auto label = fmt::format("{}", EffectType::from_genetic(gs.type));
        check(gs.coeffs.size(), design->X.cols(), label + ".coeffs");
        check(gs.u.size(), design->X.rows(), label + ".u");
    }

    check(
        state.residual().y_adj.size(),
        model.num_individuals(),
        "residual.y_adj");
}

template <typename ChainFactory>
auto Solver<ChainFactory>::resume(
    const BayesModel& model,
    bayes::BayesPrior prior,
    const std::filesystem::path& checkpoint_path,
    const MCMCObserver& observer) -> mcmc::Result
{
    mcmc::State state{model, prior};
    auto rng = read_checkpoint(checkpoint_path, state);
    validate_checkpoint(state, model);

    mcmc::Samples samples(
        model, prior, state, sample_prefix_, params_.n_records());
    notify(observer, FitPriorSetEvent{.prior = &prior});

    const infra::detail::EigenThreadGuard guard;
    omp_set_num_threads(1);

    mcmc::Context ctx{
        .model = model, .prior = prior, .state = state, .rng = rng};
    auto chain = make_chain_(ctx);
    run_impl(chain, samples, state, rng, observer);
    return finalize(std::move(samples), model, observer);
}

template <typename ChainFactory>
auto Solver<ChainFactory>::finalize(
    mcmc::Samples samples,
    const BayesModel& model,
    const MCMCObserver& observer) -> mcmc::Result
{
    samples.finalize();

    notify(
        observer,
        MCMCProgressEvent{
            .current = static_cast<size_t>(params_.n_iters),
            .total = static_cast<size_t>(params_.n_iters),
            .done = true,
        });

    mcmc::Result result(std::move(samples), model);
    result.compute();
    notify(observer, MCMCCompleteEvent{&result, &model, params_.n_records()});
    return result;
}

template <typename ChainFactory>
template <typename Chain>
auto Solver<ChainFactory>::run_impl(
    Chain& chain,
    mcmc::Samples& samples,
    mcmc::State& state,
    std::mt19937_64& rng,
    const MCMCObserver& observer) -> void
{
    for (Eigen::Index iter = 0; iter < params_.n_iters; ++iter)
    {
        chain.step();
        state.compute_heritability();

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
            samples.store(state);
        }

        const auto step = params_.checkpoint_step == 0
                              ? params_.n_iters
                              : params_.checkpoint_step;
        if (checkpoint_prefix_
            && ((iter + 1) % step == 0 || iter == params_.n_iters - 1))
        {
            write_checkpoint(state, rng, *checkpoint_prefix_);
            notify(observer, FitCheckpointSavedEvent{});
        }
    }
}

}  // namespace mcmc

}  // namespace gelex

#endif  // GELEX_ALGO_INFER_MCMC_SOLVER_H_
