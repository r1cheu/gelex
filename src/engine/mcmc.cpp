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
#include <cmath>
#include <numeric>
#include <string>
#include <utility>
#include <vector>

#include <fmt/format.h>

#include <Eigen/Core>

#include "gelex/algo/infer/mcmc/recipes.h"
#include "gelex/algo/infer/mcmc/solver.h"
#include "gelex/exception.h"
#include "gelex/infra/logging/notify.h"
#include "gelex/io/mcmc/checkpoint_reader.h"
#include "gelex/io/mcmc/result_writer.h"
#include "gelex/model/bayes/builder.h"
#include "gelex/model/bayes/model.h"
#include "gelex/model/bayes/prior.h"
#include "gelex/types/genetic_effect_type.h"

namespace gelex
{

namespace
{

using TraitRunner = mcmc::Result (*)(
    BayesModel&,
    bayes::BayesMethod,
    const mcmc::FitEngine::Config&,
    const FitObserver&);

template <auto MakeChain>
auto run_recipe(
    BayesModel& model,
    bayes::BayesMethod method,
    const mcmc::FitEngine::Config& config,
    const FitObserver& observer) -> mcmc::Result
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

auto find_spec(bayes::GeneticPrior& prior, GeneticMode mode)
    -> bayes::GeneticSpec*
{
    if (auto* s = std::get_if<bayes::GeneticSpec>(&prior.spec);
        s != nullptr && s->mode == mode)
    {
        return s;
    }
    if (auto* j = std::get_if<bayes::JointSpec>(&prior.spec); j != nullptr)
    {
        if (mode == GeneticMode::A)
        {
            return &j->additive;
        }
        if (mode == GeneticMode::D)
        {
            return &j->dominance;
        }
    }
    return nullptr;
}

auto as_eigen(const std::vector<double>& v) -> Eigen::VectorXd
{
    return Eigen::Map<const Eigen::VectorXd>(
        v.data(), static_cast<Eigen::Index>(v.size()));
}

auto apply_proportions(
    bayes::BayesMethod& method,
    GeneticMode mode,
    const std::optional<std::vector<double>>& v) -> void
{
    if (!v)
    {
        return;
    }
    for (auto& prior : method.genetics)
    {
        if (find_spec(prior, mode) != nullptr && prior.mixture)
        {
            prior.mixture->proportions.init = as_eigen(*v);
        }
    }
}

auto apply_multipliers(
    bayes::BayesMethod& method,
    GeneticMode mode,
    const std::optional<std::vector<double>>& v) -> void
{
    if (!v)
    {
        return;
    }
    for (auto& prior : method.genetics)
    {
        if (find_spec(prior, mode) == nullptr || !prior.mixture)
        {
            continue;
        }
        if (auto* sm
            = std::get_if<bayes::ScaledMixture>(&prior.mixture->strategy))
        {
            sm->multiplier = as_eigen(*v);
        }
    }
}

auto apply_dominance_sign(bayes::BayesMethod& method, double positive_prob)
    -> void
{
    for (auto& prior : method.genetics)
    {
        if (auto* spec = find_spec(prior, GeneticMode::D);
            spec != nullptr && spec->sign)
        {
            spec->sign->init
                = Eigen::VectorXd{{positive_prob, 1.0 - positive_prob}};
        }
    }
}

auto apply_prior_overrides(
    bayes::BayesMethod& method,
    const mcmc::FitEngine::Config& config) -> void
{
    apply_proportions(method, GeneticMode::A, config.pi);
    apply_proportions(method, GeneticMode::D, config.dpi);
    apply_multipliers(method, GeneticMode::A, config.multiplier);
    apply_multipliers(method, GeneticMode::D, config.dmultiplier);
    if (config.method.dominance == bayes::DominancePolicy::asymmetric)
    {
        apply_dominance_sign(method, config.positive_prob);
    }
}

}  // namespace

namespace mcmc
{

FitEngine::FitEngine(Config config) : config_(std::move(config)) {}

auto ConfigValidator::validate() const -> void
{
    check_method();
    check_mcmc_params();
    check_mixture_priors();
    check_positive_prob();
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

auto ConfigValidator::check_mixture_priors() const -> void
{
    auto check_pair = [](const std::optional<std::vector<double>>& probs,
                         const std::optional<std::vector<double>>& multipliers,
                         std::string_view probs_name,
                         std::string_view mult_name) -> void
    {
        if (probs)
        {
            if (std::ranges::any_of(*probs, [](double v) { return v < 0.0; }))
            {
                throw GelexException(
                    fmt::format("{} contains negative entries", probs_name));
            }
            const auto sum = std::accumulate(probs->begin(), probs->end(), 0.0);
            if (std::abs(sum - 1.0) > 1e-9)
            {
                throw GelexException(
                    fmt::format("{} must sum to 1, got {}", probs_name, sum));
            }
        }
        if (probs && multipliers && probs->size() != multipliers->size())
        {
            throw GelexException(
                fmt::format(
                    "{} size ({}) must match {} size ({})",
                    probs_name,
                    probs->size(),
                    mult_name,
                    multipliers->size()));
        }
    };

    check_pair(config_.pi, config_.multiplier, "pi", "multiplier");
    check_pair(config_.dpi, config_.dmultiplier, "dpi", "dmultiplier");
}

auto ConfigValidator::check_positive_prob() const -> void
{
    if (config_.positive_prob < 0.0 || config_.positive_prob > 1.0)
    {
        throw GelexException(
            fmt::format(
                "positive_prob must be in [0, 1], got {}",
                config_.positive_prob));
    }
}

auto FitEngine::run(
    PhenoPipe&& pheno,
    GenoPipe&& geno,
    const FitObserver& observer) -> void
{
    auto model = build_bayes_model(std::move(pheno), std::move(geno));
    const auto stats = compute_genetic_stats(model, config_.method);
    auto method = bayes::build_bayes_method(
        config_.method, stats, model.phenotype_variance());
    apply_prior_overrides(method, config_);
    auto runner = find_runner(config_.method);
    auto result = runner(model, std::move(method), config_, observer);

    ResultWriter writer(result, config_.bfile_prefix + ".bim");
    writer.save(config_.out_prefix);

    notify(observer, FitResultsSavedEvent{.out_prefix = config_.out_prefix});
}

}  // namespace mcmc

}  // namespace gelex
