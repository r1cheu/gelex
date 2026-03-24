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
#include <string_view>
#include <utility>
#include <vector>

#include <Eigen/Core>

#include <fmt/format.h>

#include "gelex/algo/infer/mcmc.h"
#include "gelex/exception.h"
#include "gelex/infra/logging/notify.h"
#include "gelex/infra/stats/descriptive.h"
#include "gelex/model/bayes/effects.h"
#include "gelex/model/bayes/model.h"
#include "gelex/model/bayes/prior.h"
#include "gelex/model/bayes/prior_config.h"
#include "gelex/model/bayes/trait_model.h"
#include "gelex/model/bayes/writer/result_writer.h"
#include "gelex/pipeline/geno_pipe.h"
#include "gelex/pipeline/pheno_pipe.h"

namespace gelex
{

auto FitEngine::build_model(
    PhenoPipe&& pheno,
    GenoPipe&& geno,
    const Config& config) -> BayesModel
{
    auto phenotype = std::move(pheno).take_phenotype();
    auto fixed_effects = std::move(pheno).take_fixed_effects();

    std::vector<bayes::GeneticEffect> genetics;
    genetics.emplace_back(
        GeneticKind::Add, std::move(geno).take_additive_matrix());
    if (geno.has_dominance_matrix())
    {
        genetics.emplace_back(
            GeneticKind::Dom, std::move(geno).take_dominance_matrix());
    }

    PriorSetConfig pc(config.method, detail::var(phenotype)(0));
    if (config.pi)
    {
        pc.override_proportion(GeneticKind::Add, *config.pi);
    }
    if (config.dpi)
    {
        pc.override_proportion(GeneticKind::Dom, *config.dpi);
    }
    if (config.multiplier)
    {
        pc.override_multiplier(GeneticKind::Add, *config.multiplier);
    }
    if (config.dmultiplier)
    {
        pc.override_multiplier(GeneticKind::Dom, *config.dmultiplier);
    }
    pc.override_positive_prob(config.positive_prob);
    auto priors = bayes::PriorSet::build(pc, genetics, 0);

    return BayesModel(
        std::move(phenotype),
        std::move(fixed_effects),
        std::move(genetics),
        std::move(priors));
}

namespace
{

using TraitRunner = MCMCResult (*)(
    BayesModel&,
    const MCMCParams&,
    Eigen::Index,
    std::string_view,
    const FitObserver&);

template <typename TM>
auto run_trait_model(
    BayesModel& model,
    const MCMCParams& params,
    Eigen::Index seed,
    std::string_view out_prefix,
    const FitObserver& observer) -> MCMCResult
{
    MCMC mcmc(params, TM{});
    return mcmc.run(model, seed, out_prefix, observer);
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
        throw ArgumentValidationException(
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
    auto model = build_model(std::move(pheno), std::move(geno), config_);

    notify(observer, FitModelReadyEvent{.model = &model});

    auto runner = find_runner(config_.method);
    auto result = runner(
        model, config_.mcmc_params, config_.seed, config_.out_prefix, observer);

    MCMCResultWriter writer(result, config_.bfile_prefix + ".bim");
    writer.save(config_.out_prefix);

    notify(observer, FitResultsSavedEvent{.out_prefix = config_.out_prefix});
}

}  // namespace gelex
