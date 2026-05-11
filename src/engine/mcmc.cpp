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
#include <span>
#include <string>
#include <utility>

#include <fmt/format.h>

#include "gelex/algo/infer/mcmc/recipes.h"
#include "gelex/algo/infer/mcmc/solver.h"
#include "gelex/exception.h"
#include "gelex/infra/logging/notify.h"
#include "gelex/io/mcmc/checkpoint_reader.h"
#include "gelex/io/mcmc/result_writer.h"
#include "gelex/model/bayes/algorithm_shape.h"
#include "gelex/model/bayes/bayes_policy.h"
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

struct DispatchKey
{
    BayesBase base{};
    bayes::AlgorithmShape shape{};
    bayes::DominancePolicy dominance{};
    bool estimate_pi{};

    constexpr auto operator==(const DispatchKey&) const -> bool = default;
};

// clang-format off
constexpr auto kSym = bayes::DominancePolicy::symmetric;
using Shape = bayes::AlgorithmShape;

constexpr auto kTraitRunners = std::array<
    std::pair<DispatchKey, TraitRunner>, 21>{{
    {{BayesBase::A,  Shape::a_only,         kSym, false}, &run_recipe<mcmc::make_bayes_a_chain<Shape::a_only>>},
    {{BayesBase::A,  Shape::d_only,         kSym, false}, &run_recipe<mcmc::make_bayes_a_chain<Shape::d_only>>},
    {{BayesBase::A,  Shape::ad_independent, kSym, false}, &run_recipe<mcmc::make_bayes_a_chain<Shape::ad_independent>>},
    {{BayesBase::B,  Shape::a_only,         kSym, false}, &run_recipe<mcmc::make_bayes_b_chain<Shape::a_only>>},
    {{BayesBase::B,  Shape::d_only,         kSym, false}, &run_recipe<mcmc::make_bayes_b_chain<Shape::d_only>>},
    {{BayesBase::B,  Shape::ad_independent, kSym, false}, &run_recipe<mcmc::make_bayes_b_chain<Shape::ad_independent>>},
    {{BayesBase::B,  Shape::a_only,         kSym, true},  &run_recipe<mcmc::make_bayes_bpi_chain<Shape::a_only>>},
    {{BayesBase::B,  Shape::d_only,         kSym, true},  &run_recipe<mcmc::make_bayes_bpi_chain<Shape::d_only>>},
    {{BayesBase::B,  Shape::ad_independent, kSym, true},  &run_recipe<mcmc::make_bayes_bpi_chain<Shape::ad_independent>>},
    {{BayesBase::C,  Shape::a_only,         kSym, false}, &run_recipe<mcmc::make_bayes_c_chain<Shape::a_only>>},
    {{BayesBase::C,  Shape::d_only,         kSym, false}, &run_recipe<mcmc::make_bayes_c_chain<Shape::d_only>>},
    {{BayesBase::C,  Shape::ad_independent, kSym, false}, &run_recipe<mcmc::make_bayes_c_chain<Shape::ad_independent>>},
    {{BayesBase::C,  Shape::a_only,         kSym, true},  &run_recipe<mcmc::make_bayes_cpi_chain<Shape::a_only>>},
    {{BayesBase::C,  Shape::d_only,         kSym, true},  &run_recipe<mcmc::make_bayes_cpi_chain<Shape::d_only>>},
    {{BayesBase::C,  Shape::ad_independent, kSym, true},  &run_recipe<mcmc::make_bayes_cpi_chain<Shape::ad_independent>>},
    {{BayesBase::R,  Shape::a_only,         kSym, false}, &run_recipe<mcmc::make_bayes_r_chain<Shape::a_only>>},
    {{BayesBase::R,  Shape::d_only,         kSym, false}, &run_recipe<mcmc::make_bayes_r_chain<Shape::d_only>>},
    {{BayesBase::R,  Shape::ad_independent, kSym, false}, &run_recipe<mcmc::make_bayes_r_chain<Shape::ad_independent>>},
    {{BayesBase::RR, Shape::a_only,         kSym, false}, &run_recipe<mcmc::make_bayes_rr_chain<Shape::a_only>>},
    {{BayesBase::RR, Shape::d_only,         kSym, false}, &run_recipe<mcmc::make_bayes_rr_chain<Shape::d_only>>},
    {{BayesBase::RR, Shape::ad_independent, kSym, false}, &run_recipe<mcmc::make_bayes_rr_chain<Shape::ad_independent>>},
}};
// clang-format on

auto find_runner(
    const bayes::BayesConfig& method,
    std::span<const GeneticMode> requested) -> TraitRunner
{
    const auto shape
        = bayes::resolve_shape(bayes::policy_for(method.base), requested);
    const DispatchKey key{
        method.base, shape, method.dominance, method.estimate_pi};

    const auto* it = std::ranges::find(
        kTraitRunners, key, &std::pair<DispatchKey, TraitRunner>::first);
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
    if (!bayes::is_valid_method(config_.method, config_.requested_effects))
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
    auto runner = find_runner(config_.method, config_.requested_effects);
    auto result = runner(model, std::move(method), config_, observer);

    ResultWriter writer(result, config_.bfile_prefix + ".bim");
    writer.save(config_.out_prefix);

    notify(observer, FitResultsSavedEvent{.out_prefix = config_.out_prefix});
}

}  // namespace mcmc

}  // namespace gelex
