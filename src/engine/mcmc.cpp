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

#include <algorithm>
#include <array>
#include <string>
#include <utility>

#include <fmt/format.h>

#include "gelex/algo/infer/mcmc/recipes.h"
#include "gelex/algo/infer/mcmc/solver.h"
#include "gelex/exception.h"
#include "gelex/infra/logging/notify.h"
#include "gelex/io/mcmc/checkpoint_reader.h"
#include "gelex/io/mcmc/result_writer.h"
#include "gelex/model/bayes/model.h"
#include "gelex/types/genetic_effect_type.h"

namespace gelex
{

namespace
{

using TraitRunner = mcmc::Result (*)(
    const BayesModel&,
    bayes::BayesMethod,
    const mcmc::Engine::Config&,
    const MCMCObserver&);

template <auto MakeChain>
auto run_recipe(
    const BayesModel& model,
    bayes::BayesMethod method,
    const mcmc::Engine::Config& config,
    const MCMCObserver& observer) -> mcmc::Result
{
    mcmc::Solver mcmc(
        config.mcmc_params,
        MakeChain,
        std::string(config.out_prefix),
        std::string(config.out_prefix));
    if (config.resume_path)
    {
        auto checkpoint = read_checkpoint(*config.resume_path);
        return mcmc.resume(model, std::move(checkpoint), observer);
    }
    return mcmc.run(model, std::move(method), config.seed, observer);
}

// clang-format off
constexpr auto kSym = bayes::DominancePolicy::symmetric;

constexpr auto kTraitRunners = std::array<
    std::pair<bayes::BayesConfig, TraitRunner>, 21>{{
    {{BayesBase::A,  GeneticMode::A,  kSym, false}, &run_recipe<mcmc::make_bayes_a_chain<GeneticMode::A>>},
    {{BayesBase::A,  GeneticMode::D,  kSym, false}, &run_recipe<mcmc::make_bayes_a_chain<GeneticMode::D>>},
    {{BayesBase::A,  GeneticMode::AD, kSym, false}, &run_recipe<mcmc::make_bayes_a_chain<GeneticMode::AD>>},
    {{BayesBase::B,  GeneticMode::A,  kSym, false}, &run_recipe<mcmc::make_bayes_b_chain<GeneticMode::A>>},
    {{BayesBase::B,  GeneticMode::D,  kSym, false}, &run_recipe<mcmc::make_bayes_b_chain<GeneticMode::D>>},
    {{BayesBase::B,  GeneticMode::AD, kSym, false}, &run_recipe<mcmc::make_bayes_b_chain<GeneticMode::AD>>},
    {{BayesBase::B,  GeneticMode::A,  kSym, true},  &run_recipe<mcmc::make_bayes_bpi_chain<GeneticMode::A>>},
    {{BayesBase::B,  GeneticMode::D,  kSym, true},  &run_recipe<mcmc::make_bayes_bpi_chain<GeneticMode::D>>},
    {{BayesBase::B,  GeneticMode::AD, kSym, true},  &run_recipe<mcmc::make_bayes_bpi_chain<GeneticMode::AD>>},
    {{BayesBase::C,  GeneticMode::A,  kSym, false}, &run_recipe<mcmc::make_bayes_c_chain<GeneticMode::A>>},
    {{BayesBase::C,  GeneticMode::D,  kSym, false}, &run_recipe<mcmc::make_bayes_c_chain<GeneticMode::D>>},
    {{BayesBase::C,  GeneticMode::AD, kSym, false}, &run_recipe<mcmc::make_bayes_c_chain<GeneticMode::AD>>},
    {{BayesBase::C,  GeneticMode::A,  kSym, true},  &run_recipe<mcmc::make_bayes_cpi_chain<GeneticMode::A>>},
    {{BayesBase::C,  GeneticMode::D,  kSym, true},  &run_recipe<mcmc::make_bayes_cpi_chain<GeneticMode::D>>},
    {{BayesBase::C,  GeneticMode::AD, kSym, true},  &run_recipe<mcmc::make_bayes_cpi_chain<GeneticMode::AD>>},
    {{BayesBase::R,  GeneticMode::A,  kSym, false}, &run_recipe<mcmc::make_bayes_r_chain<GeneticMode::A>>},
    {{BayesBase::R,  GeneticMode::D,  kSym, false}, &run_recipe<mcmc::make_bayes_r_chain<GeneticMode::D>>},
    {{BayesBase::R,  GeneticMode::AD, kSym, false}, &run_recipe<mcmc::make_bayes_r_chain<GeneticMode::AD>>},
    {{BayesBase::RR, GeneticMode::A,  kSym, false}, &run_recipe<mcmc::make_bayes_rr_chain<GeneticMode::A>>},
    {{BayesBase::RR, GeneticMode::D,  kSym, false}, &run_recipe<mcmc::make_bayes_rr_chain<GeneticMode::D>>},
    {{BayesBase::RR, GeneticMode::AD, kSym, false}, &run_recipe<mcmc::make_bayes_rr_chain<GeneticMode::AD>>},
}};
// clang-format on

auto find_runner(bayes::BayesConfig method) -> TraitRunner
{
    const auto* it = std::ranges::find(
        kTraitRunners,
        method,
        &std::pair<bayes::BayesConfig, TraitRunner>::first);
    if (it == kTraitRunners.end())
    {
        throw GelexException(
            fmt::format("Unsupported Bayes method: {}", method));
    }
    return it->second;
}

}  // namespace

namespace mcmc
{

Engine::Engine(Config config) : config_(std::move(config)) {}

auto ConfigValidator::validate() const -> void
{
    check_method();
    check_mcmc_params();
}

auto ConfigValidator::check_method() const -> void
{
    if (!bayes::is_valid_method(config_.method))
    {
        throw GelexException(
            fmt::format("invalid method combination: {}", config_.method));
    }
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
    bayes::BayesMethod method,
    const MCMCObserver& observer) -> void
{
    auto runner = find_runner(config_.method);
    auto result = runner(model, std::move(method), config_, observer);

    ResultWriter writer(result, config_.bfile_prefix + ".bim");
    writer.save(config_.out_prefix);

    notify(observer, FitResultsSavedEvent{.out_prefix = config_.out_prefix});
}

}  // namespace mcmc

}  // namespace gelex
