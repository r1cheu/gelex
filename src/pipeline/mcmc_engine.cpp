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

#include "gelex/pipeline/mcmc_engine.h"

#include <algorithm>
#include <array>
#include <string>
#include <utility>

#include <fmt/format.h>

#include "gelex/algo/infer/mcmc/solver.h"
#include "gelex/algo/infer/mcmc/trait_model.h"
#include "gelex/exception.h"
#include "gelex/infra/logging/notify.h"
#include "gelex/model/bayes/model.h"
#include "gelex/model/bayes/reader/checkpoint_reader.h"
#include "gelex/model/bayes/writer/mcmc_result_writer.h"
#include "gelex/pipeline/model_builder.h"

namespace gelex
{

namespace
{

using TraitRunner = mcmc::Result (*)(
    BayesModel&,
    const bayes::Priors&,
    const mcmc::FitEngine::Config&,
    const FitObserver&);

template <typename TM>
auto run_trait_model(
    BayesModel& model,
    const bayes::Priors& priors,
    const mcmc::FitEngine::Config& config,
    const FitObserver& observer) -> mcmc::Result
{
    mcmc::Solver mcmc(config.mcmc_params, TM{}, std::string(config.out_prefix));
    if (config.resume_path)
    {
        auto checkpoint = read_checkpoint(*config.resume_path);
        return mcmc.resume(
            model, std::move(checkpoint), config.out_prefix, observer);
    }
    return mcmc.run(model, priors, config.seed, config.out_prefix, observer);
}

// clang-format off
constexpr auto kTraitRunners = std::array<
    std::pair<BayesMethodConfig, TraitRunner>, 15>{{
    {{BayesBase::A,  false, false, false}, &run_trait_model<mcmc::A>},
    {{BayesBase::A,  true,  false, false}, &run_trait_model<mcmc::Ad>},
    {{BayesBase::B,  false, false, false}, &run_trait_model<mcmc::B>},
    {{BayesBase::B,  false, false, true},  &run_trait_model<mcmc::Bpi>},
    {{BayesBase::B,  true,  false, false}, &run_trait_model<mcmc::Bd>},
    {{BayesBase::B,  true,  false, true},  &run_trait_model<mcmc::Bdpi>},
    {{BayesBase::C,  false, false, false}, &run_trait_model<mcmc::C>},
    {{BayesBase::C,  false, false, true},  &run_trait_model<mcmc::Cpi>},
    {{BayesBase::C,  true,  false, false}, &run_trait_model<mcmc::Cd>},
    {{BayesBase::C,  true,  false, true},  &run_trait_model<mcmc::Cdpi>},
    {{BayesBase::R,  false, false, false}, &run_trait_model<mcmc::R>},
    {{BayesBase::R,  true,  false, false}, &run_trait_model<mcmc::Rd>},
    {{BayesBase::R,  true,  true,  false}, &run_trait_model<mcmc::RdAt>},
    {{BayesBase::RR, false, false, false}, &run_trait_model<mcmc::RR>},
    {{BayesBase::RR, true,  false, false}, &run_trait_model<mcmc::RRd>},
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
