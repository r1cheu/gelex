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

#include "gelex/pipeline/fit_engine.h"

#include <algorithm>
#include <array>
#include <string>
#include <utility>
#include <vector>

#include <fmt/format.h>

#include "gelex/algo/infer/mcmc.h"
#include "gelex/exception.h"
#include "gelex/infra/logging/notify.h"
#include "gelex/model/bayes/effects.h"
#include "gelex/model/bayes/model.h"
#include "gelex/model/bayes/prior.h"
#include "gelex/model/bayes/prior_config.h"
#include "gelex/model/bayes/reader/checkpoint_reader.h"
#include "gelex/model/bayes/trait_model.h"
#include "gelex/model/bayes/writer/result_writer.h"
#include "gelex/pipeline/geno_pipe.h"
#include "gelex/pipeline/pheno_pipe.h"

namespace gelex
{

auto FitEngine::build_model(PhenoPipe&& pheno, GenoPipe&& geno) -> BayesModel
{
    auto phenotype = std::move(pheno).take_phenotype();
    auto fixed_effects = std::move(pheno).take_fixed_effects();

    std::vector<bayes::GeneticEffect> genetics;
    genetics.emplace_back(
        GeneticMode::A, std::move(geno).take_additive_matrix());
    if (geno.has_dominance_matrix())
    {
        genetics.emplace_back(
            GeneticMode::D, std::move(geno).take_dominance_matrix());
    }

    return BayesModel(
        std::move(phenotype), std::move(fixed_effects), std::move(genetics));
}

auto FitEngine::build_priors(const Config& config, const BayesModel& model)
    -> bayes::Priors
{
    PriorSetConfig pc(config.method, model.phenotype_variance());
    if (config.pi)
    {
        pc.override_proportion(GeneticMode::A, *config.pi);
    }
    if (config.dpi)
    {
        pc.override_proportion(GeneticMode::D, *config.dpi);
    }
    if (config.multiplier)
    {
        pc.override_multiplier(GeneticMode::A, *config.multiplier);
    }
    if (config.dmultiplier)
    {
        pc.override_multiplier(GeneticMode::D, *config.dmultiplier);
    }
    pc.override_positive_prob(config.positive_prob);
    return bayes::Priors(pc, model.genetics(), 0);
}

namespace
{

using TraitRunner = MCMCResult (*)(
    BayesModel&,
    const bayes::Priors&,
    const FitEngine::Config&,
    const FitObserver&);

template <typename TM>
auto run_trait_model(
    BayesModel& model,
    const bayes::Priors& priors,
    const FitEngine::Config& config,
    const FitObserver& observer) -> MCMCResult
{
    MCMC mcmc(config.mcmc_params, TM{}, std::string(config.out_prefix));
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
    {{BayesBase::A,  false, false, false}, &run_trait_model<BayesA>},
    {{BayesBase::A,  true,  false, false}, &run_trait_model<BayesAd>},
    {{BayesBase::B,  false, false, false}, &run_trait_model<BayesB>},
    {{BayesBase::B,  false, false, true},  &run_trait_model<BayesBpi>},
    {{BayesBase::B,  true,  false, false}, &run_trait_model<BayesBd>},
    {{BayesBase::B,  true,  false, true},  &run_trait_model<BayesBdpi>},
    {{BayesBase::C,  false, false, false}, &run_trait_model<BayesC>},
    {{BayesBase::C,  false, false, true},  &run_trait_model<BayesCpi>},
    {{BayesBase::C,  true,  false, false}, &run_trait_model<BayesCd>},
    {{BayesBase::C,  true,  false, true},  &run_trait_model<BayesCdpi>},
    {{BayesBase::R,  false, false, false}, &run_trait_model<BayesR>},
    {{BayesBase::R,  true,  false, false}, &run_trait_model<BayesRd>},
    {{BayesBase::R,  true,  true,  false}, &run_trait_model<BayesRdAt>},
    {{BayesBase::RR, false, false, false}, &run_trait_model<BayesRR>},
    {{BayesBase::RR, true,  false, false}, &run_trait_model<BayesRRd>},
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

FitEngine::FitEngine(Config config) : config_(std::move(config)) {}

auto FitEngine::run(
    PhenoPipe&& pheno,
    GenoPipe&& geno,
    const FitObserver& observer) -> void
{
    auto model = build_model(std::move(pheno), std::move(geno));
    auto priors = build_priors(config_, model);
    auto runner = find_runner(config_.method);
    auto result = runner(model, priors, config_, observer);

    MCMCResultWriter writer(result, config_.bfile_prefix + ".bim");
    writer.save(config_.out_prefix);

    notify(observer, FitResultsSavedEvent{.out_prefix = config_.out_prefix});
}

}  // namespace gelex
