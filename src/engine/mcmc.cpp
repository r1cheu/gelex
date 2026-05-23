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
#include <span>
#include <string>
#include <utility>

#include <fmt/format.h>

#include "gelex/algo/infer/mcmc/recipes.h"
#include "gelex/algo/infer/mcmc/solver.h"
#include "gelex/exception.h"
#include "gelex/infra/logging/notify.h"
#include "gelex/io/mcmc/result_writer.h"
#include "gelex/model/bayes/model.h"
#include "gelex/model/bayes/prior.h"
#include "gelex/model/bayes/recipe.h"
#include "gelex/types/genetic_effect_type.h"

namespace gelex
{

namespace
{

using TraitRunner = mcmc::Result (*)(
    const BayesModel&,
    bayes::BayesPrior,
    const mcmc::Engine::Config&,
    const MCMCObserver&);

template <auto MakeChain>
auto run_recipe(
    const BayesModel& model,
    bayes::BayesPrior prior,
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
        return mcmc.resume(
            model, std::move(prior), *config.resume_path, observer);
    }
    return mcmc.run(model, std::move(prior), config.seed, observer);
}

auto resolve_shape(std::span<const GeneticMode> modes) -> mcmc::GeneticShape
{
    if (modes.size() == 1 && modes[0] == GeneticMode::A)
    {
        return mcmc::GeneticShape::a_only;
    }
    if (modes.size() == 1 && modes[0] == GeneticMode::D)
    {
        return mcmc::GeneticShape::d_only;
    }
    if (modes.size() == 2 && modes[0] == GeneticMode::A
        && modes[1] == GeneticMode::D)
    {
        return mcmc::GeneticShape::ad_independent;
    }
    throw GelexException("MCMC engine supports modes {A}, {D}, or {A, D}");
}

template <bayes::BayesRecipePreset Preset, mcmc::GeneticShape Shape>
consteval auto runner_for() -> TraitRunner
{
    if constexpr (Preset == bayes::BayesRecipePreset::RR)
    {
        return &run_recipe<mcmc::make_bayes_rr_chain<Shape>>;
    }
    else if constexpr (Preset == bayes::BayesRecipePreset::A)
    {
        return &run_recipe<mcmc::make_bayes_a_chain<Shape>>;
    }
    else if constexpr (Preset == bayes::BayesRecipePreset::B)
    {
        return &run_recipe<mcmc::make_bayes_bpi_chain<Shape>>;
    }
    else if constexpr (Preset == bayes::BayesRecipePreset::C)
    {
        return &run_recipe<mcmc::make_bayes_cpi_chain<Shape>>;
    }
    else if constexpr (Preset == bayes::BayesRecipePreset::R)
    {
        return &run_recipe<mcmc::make_bayes_r_chain<Shape>>;
    }
}

template <bayes::BayesRecipePreset Preset>
auto runner_for_shape(mcmc::GeneticShape shape) -> TraitRunner
{
    switch (shape)
    {
        case mcmc::GeneticShape::a_only:
            return runner_for<Preset, mcmc::GeneticShape::a_only>();
        case mcmc::GeneticShape::d_only:
            return runner_for<Preset, mcmc::GeneticShape::d_only>();
        case mcmc::GeneticShape::ad_independent:
            return runner_for<Preset, mcmc::GeneticShape::ad_independent>();
    }
    std::unreachable();
}

auto find_runner(
    bayes::BayesRecipePreset preset,
    std::span<const GeneticMode> modes) -> TraitRunner
{
    if (preset == bayes::BayesRecipePreset::CD)
    {
        throw GelexException(
            "CD MCMC runtime is not implemented for BayesPrior yet");
    }

    const auto shape = resolve_shape(modes);
    switch (preset)
    {
        case bayes::BayesRecipePreset::RR:
            return runner_for_shape<bayes::BayesRecipePreset::RR>(shape);
        case bayes::BayesRecipePreset::A:
            return runner_for_shape<bayes::BayesRecipePreset::A>(shape);
        case bayes::BayesRecipePreset::B:
            return runner_for_shape<bayes::BayesRecipePreset::B>(shape);
        case bayes::BayesRecipePreset::C:
            return runner_for_shape<bayes::BayesRecipePreset::C>(shape);
        case bayes::BayesRecipePreset::R:
            return runner_for_shape<bayes::BayesRecipePreset::R>(shape);
        case bayes::BayesRecipePreset::CD:
            break;
    }
    std::unreachable();
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
    if (config_.preset == bayes::BayesRecipePreset::CD)
    {
        throw GelexException(
            "CD MCMC runtime is not implemented for BayesPrior yet");
    }
    static_cast<void>(resolve_shape(config_.recipe_config.modes));
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
    auto runner = find_runner(config_.preset, config_.recipe_config.modes);
    auto result = runner(model, std::move(prior), config_, observer);

    ResultWriter writer(result, config_.bfile_prefix + ".bim");
    writer.save(config_.out_prefix);

    notify(observer, FitResultsSavedEvent{.out_prefix = config_.out_prefix});
}

}  // namespace mcmc

}  // namespace gelex
