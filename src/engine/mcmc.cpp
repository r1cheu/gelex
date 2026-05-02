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
#include "gelex/model/bayes/builder.h"
#include "gelex/model/bayes/model.h"
#include "gelex/types/genetic_effect_type.h"

namespace gelex
{

namespace
{

using TraitRunner = mcmc::Result (*)(
    BayesModel&,
    const bayes::Priors&,
    const mcmc::FitEngine::Config&,
    const FitObserver&);

template <auto MakeChain>
auto run_recipe(
    BayesModel& model,
    const bayes::Priors& priors,
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
    return mcmc.run(model, priors, config.seed, observer);
}

// clang-format off
constexpr auto kTraitRunners = std::array<
    std::pair<BayesMethodConfig, TraitRunner>, 21>{{
    {{BayesBase::A,  GeneticMode::A,  false, false}, &run_recipe<mcmc::make_bayes_a_chain<GeneticMode::A>>},
    {{BayesBase::A,  GeneticMode::D,  false, false}, &run_recipe<mcmc::make_bayes_a_chain<GeneticMode::D>>},
    {{BayesBase::A,  GeneticMode::AD, false, false}, &run_recipe<mcmc::make_bayes_a_chain<GeneticMode::AD>>},
    {{BayesBase::B,  GeneticMode::A,  false, false}, &run_recipe<mcmc::make_bayes_b_chain<GeneticMode::A>>},
    {{BayesBase::B,  GeneticMode::D,  false, false}, &run_recipe<mcmc::make_bayes_b_chain<GeneticMode::D>>},
    {{BayesBase::B,  GeneticMode::AD, false, false}, &run_recipe<mcmc::make_bayes_b_chain<GeneticMode::AD>>},
    {{BayesBase::B,  GeneticMode::A,  false, true},  &run_recipe<mcmc::make_bayes_bpi_chain<GeneticMode::A>>},
    {{BayesBase::B,  GeneticMode::D,  false, true},  &run_recipe<mcmc::make_bayes_bpi_chain<GeneticMode::D>>},
    {{BayesBase::B,  GeneticMode::AD, false, true},  &run_recipe<mcmc::make_bayes_bpi_chain<GeneticMode::AD>>},
    {{BayesBase::C,  GeneticMode::A,  false, false}, &run_recipe<mcmc::make_bayes_c_chain<GeneticMode::A>>},
    {{BayesBase::C,  GeneticMode::D,  false, false}, &run_recipe<mcmc::make_bayes_c_chain<GeneticMode::D>>},
    {{BayesBase::C,  GeneticMode::AD, false, false}, &run_recipe<mcmc::make_bayes_c_chain<GeneticMode::AD>>},
    {{BayesBase::C,  GeneticMode::A,  false, true},  &run_recipe<mcmc::make_bayes_cpi_chain<GeneticMode::A>>},
    {{BayesBase::C,  GeneticMode::D,  false, true},  &run_recipe<mcmc::make_bayes_cpi_chain<GeneticMode::D>>},
    {{BayesBase::C,  GeneticMode::AD, false, true},  &run_recipe<mcmc::make_bayes_cpi_chain<GeneticMode::AD>>},
    {{BayesBase::R,  GeneticMode::A,  false, false}, &run_recipe<mcmc::make_bayes_r_chain<GeneticMode::A>>},
    {{BayesBase::R,  GeneticMode::D,  false, false}, &run_recipe<mcmc::make_bayes_r_chain<GeneticMode::D>>},
    {{BayesBase::R,  GeneticMode::AD, false, false}, &run_recipe<mcmc::make_bayes_r_chain<GeneticMode::AD>>},
    {{BayesBase::RR, GeneticMode::A,  false, false}, &run_recipe<mcmc::make_bayes_rr_chain<GeneticMode::A>>},
    {{BayesBase::RR, GeneticMode::D,  false, false}, &run_recipe<mcmc::make_bayes_rr_chain<GeneticMode::D>>},
    {{BayesBase::RR, GeneticMode::AD, false, false}, &run_recipe<mcmc::make_bayes_rr_chain<GeneticMode::AD>>},
}};
// clang-format on

auto find_runner(BayesMethodConfig method) -> TraitRunner
{
    const auto* it = std::ranges::find(
        kTraitRunners,
        method,
        &std::pair<BayesMethodConfig, TraitRunner>::first);
    if (it == kTraitRunners.end())
    {
        throw GelexException(
            fmt::format("Unsupported Bayes method: {}", method));
    }
    return it->second;
}

}  // namespace

mcmc::FitEngine::FitEngine(Config config) : config_(std::move(config)) {}

auto mcmc::FitEngine::run(
    PhenoPipe&& pheno,
    GenoPipe&& geno,
    const FitObserver& observer) -> void
{
    auto model = build_bayes_model(std::move(pheno), std::move(geno));
    auto priors = build_bayes_priors(
        PriorOverrides{
            .method = config_.method,
            .phenotype_variance = model.phenotype_variance(),
            .pi = config_.pi,
            .dpi = config_.dpi,
            .multiplier = config_.multiplier,
            .dmultiplier = config_.dmultiplier,
            .positive_prob = config_.positive_prob,
        },
        model);
    auto runner = find_runner(config_.method);
    auto result = runner(model, priors, config_, observer);

    mcmc::ResultWriter writer(result, config_.bfile_prefix + ".bim");
    writer.save(config_.out_prefix);

    notify(observer, FitResultsSavedEvent{.out_prefix = config_.out_prefix});
}

}  // namespace gelex
