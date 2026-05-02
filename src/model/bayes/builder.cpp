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

#include "gelex/model/bayes/builder.h"

#include <utility>
#include <vector>

#include "gelex/data/pipe/geno.h"
#include "gelex/data/pipe/pheno.h"
#include "gelex/model/bayes/effects.h"
#include "gelex/model/bayes/model.h"
#include "gelex/model/bayes/prior_config.h"
#include "gelex/types/genetic_effect_type.h"

namespace gelex
{

auto build_bayes_model(PhenoPipe&& pheno, GenoPipe&& geno) -> BayesModel
{
    auto phenotype = std::move(pheno).take_phenotype();
    auto fixed_effects = std::move(pheno).take_fixed_effects();

    std::vector<bayes::GeneticEffect> genetics;
    if (geno.has_additive_matrix())
    {
        genetics.emplace_back(
            GeneticMode::A, std::move(geno).take_additive_matrix());
    }
    if (geno.has_dominance_matrix())
    {
        genetics.emplace_back(
            GeneticMode::D, std::move(geno).take_dominance_matrix());
    }

    return BayesModel(
        std::move(phenotype), std::move(fixed_effects), std::move(genetics));
}

auto build_bayes_priors(
    const PriorOverrides& overrides,
    const BayesModel& model) -> bayes::Priors
{
    PriorSetConfig pc(overrides.method, overrides.phenotype_variance);
    if (overrides.pi)
    {
        pc.override_proportion(GeneticMode::A, *overrides.pi);
    }
    if (overrides.dpi)
    {
        pc.override_proportion(GeneticMode::D, *overrides.dpi);
    }
    if (overrides.multiplier)
    {
        pc.override_multiplier(GeneticMode::A, *overrides.multiplier);
    }
    if (overrides.dmultiplier)
    {
        pc.override_multiplier(GeneticMode::D, *overrides.dmultiplier);
    }
    pc.override_positive_prob(overrides.positive_prob);
    return bayes::Priors(pc, model.genetics(), 0);
}

}  // namespace gelex
