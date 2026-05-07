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

#include <optional>
#include <string>
#include <string_view>
#include <utility>

#include <fmt/format.h>
#include <omp.h>
#include <Eigen/Core>

#include "gelex/algo/detail/posterior_calculator.h"
#include "gelex/algo/infer/mcmc/checkpoint.h"
#include "gelex/algo/infer/mcmc/context.h"
#include "gelex/algo/infer/mcmc/result.h"
#include "gelex/algo/infer/mcmc/samples.h"
#include "gelex/algo/infer/params.h"
#include "gelex/data/genotype/storage.h"
#include "gelex/exception.h"
#include "gelex/infra/detail/eigen_thread_guard.h"
#include "gelex/infra/logging/fit_event.h"
#include "gelex/infra/logging/notify.h"
#include "gelex/io/mcmc/checkpoint_writer.h"
#include "gelex/model/bayes/method.h"
#include "gelex/model/bayes/model.h"

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
        bayes::BayesMethod method,
        Eigen::Index seed = 42,
        const FitObserver& observer = {}) -> mcmc::Result;

    auto resume(
        const BayesModel& model,
        Checkpoint checkpoint,
        const FitObserver& observer = {}) -> mcmc::Result;

   private:
    static void validate_checkpoint(
        const mcmc::State& state,
        const BayesModel& model);

    template <typename Chain>
    auto run_impl(
        Chain& chain,
        const bayes::BayesMethod& method,
        mcmc::Samples& samples,
        mcmc::State& state,
        std::mt19937_64& rng,
        const FitObserver& observer) -> void;

    auto finalize(
        mcmc::Samples samples,
        const BayesModel& model,
        const FitObserver& observer) -> mcmc::Result;

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
    bayes::BayesMethod method,
    Eigen::Index seed,
    const FitObserver& observer) -> mcmc::Result
{
    mcmc::Samples samples(model, method, sample_prefix_, params_.n_records);
    notify(observer, FitMethodSetEvent{.method = &method});

    mcmc::State state{model, method};
    std::mt19937_64 rng(seed);

    const infra::detail::EigenThreadGuard guard;
    omp_set_num_threads(1);

    mcmc::Context ctx{
        .model = model, .method = method, .state = state, .rng = rng};
    auto chain = make_chain_(ctx);
    run_impl(chain, method, samples, state, rng, observer);
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
                "checkpoint random effect count mismatch: expected {}, got {}",
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
        const auto* effect = model.genetic(gs.type);
        if (effect == nullptr)
        {
            throw GelexException(
                fmt::format(
                    "checkpoint has genetic type {} not in model",
                    EffectType::from_genetic(gs.type)));
        }
        const auto label = fmt::format("{}", EffectType::from_genetic(gs.type));
        check(gs.coeffs.size(), bayes::get_cols(effect->X), label + ".coeffs");
        check(gs.u.size(), bayes::get_rows(effect->X), label + ".u");
    }

    check(
        state.residual().y_adj.size(),
        model.num_individuals(),
        "residual.y_adj");
}

template <typename ChainFactory>
auto Solver<ChainFactory>::resume(
    const BayesModel& model,
    Checkpoint checkpoint,
    const FitObserver& observer) -> mcmc::Result
{
    validate_checkpoint(checkpoint.state, model);

    mcmc::Samples samples(
        model, checkpoint.method, sample_prefix_, params_.n_records);
    notify(observer, FitMethodSetEvent{.method = &checkpoint.method});

    const infra::detail::EigenThreadGuard guard;
    omp_set_num_threads(1);

    mcmc::Context ctx{
        .model = model,
        .method = checkpoint.method,
        .state = checkpoint.state,
        .rng = checkpoint.rng};
    auto chain = make_chain_(ctx);
    run_impl(
        chain,
        checkpoint.method,
        samples,
        checkpoint.state,
        checkpoint.rng,
        observer);
    return finalize(std::move(samples), model, observer);
}

template <typename ChainFactory>
auto Solver<ChainFactory>::finalize(
    mcmc::Samples samples,
    const BayesModel& model,
    const FitObserver& observer) -> mcmc::Result
{
    samples.finalize();

    notify(
        observer,
        FitMCMCProgressEvent{
            .current = static_cast<size_t>(params_.n_iters),
            .total = static_cast<size_t>(params_.n_iters),
            .done = true,
        });

    mcmc::Result result(std::move(samples), model);
    result.compute();
    notify(observer, FitMCMCCompleteEvent{&result, &model, params_.n_records});
    return result;
}

template <typename ChainFactory>
template <typename Chain>
auto Solver<ChainFactory>::run_impl(
    Chain& chain,
    const bayes::BayesMethod& method,
    mcmc::Samples& samples,
    mcmc::State& state,
    std::mt19937_64& rng,
    const FitObserver& observer) -> void
{
    for (Eigen::Index iter = 0; iter < params_.n_iters; ++iter)
    {
        chain.step();
        state.compute_heritability();

        notify(
            observer,
            FitMCMCProgressEvent{
                .current = static_cast<size_t>(iter + 1),
                .total = static_cast<size_t>(params_.n_iters),
                .state = &state,
            });

        if (iter >= params_.n_burn_in
            && (iter + 1 - params_.n_burn_in) % params_.n_thin == 0)
        {
            samples.store(state);
        }

        if (checkpoint_prefix_
            && ((iter + 1) % params_.checkpoint_step == 0
                || iter == params_.n_iters - 1))
        {
            write_checkpoint(state, rng, method, *checkpoint_prefix_);
            notify(observer, FitCheckpointSavedEvent{});
        }
    }
}

}  // namespace mcmc

}  // namespace gelex

#endif  // GELEX_ALGO_INFER_MCMC_SOLVER_H_
