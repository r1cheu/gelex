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

#include <algorithm>
#include <utility>
#include <vector>

#include <fmt/format.h>

#include "gelex/data/genotype/storage.h"
#include "gelex/data/pipe/geno.h"
#include "gelex/data/pipe/pheno.h"
#include "gelex/exception.h"
#include "gelex/infra/stats/descriptive.h"
#include "gelex/model/bayes/bayes_policy.h"
#include "gelex/model/bayes/effects.h"
#include "gelex/model/bayes/method.h"
#include "gelex/model/bayes/model.h"
#include "gelex/model/bayes/prior.h"
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

auto compute_genetic_stats(
    const BayesModel& model,
    const bayes::BayesConfig& config) -> std::vector<bayes::GeneticStats>
{
    constexpr double kAdditiveHeritability = 0.5;
    constexpr double kDominanceHeritability = 0.2;
    auto heritability_for = [](GeneticMode m)
    {
        return m == GeneticMode::D ? kDominanceHeritability
                                   : kAdditiveHeritability;
    };

    std::vector<GeneticMode> modes;
    switch (config.mode)
    {
        case GeneticMode::A:
            modes = {GeneticMode::A};
            break;
        case GeneticMode::D:
            modes = {GeneticMode::D};
            break;
        case GeneticMode::AD:
            modes = {GeneticMode::A, GeneticMode::D};
            break;
        case GeneticMode::kCount:
            throw GelexException("Invalid GeneticMode: kCount");
    }

    const auto& genetics = model.genetics();
    std::vector<bayes::GeneticStats> stats;
    stats.reserve(modes.size());
    for (auto mode : modes)
    {
        auto it
            = std::ranges::find(genetics, mode, &bayes::GeneticEffect::type);
        if (it == genetics.end())
        {
            throw GelexException(
                fmt::format(
                    "BayesModel missing genetic effect for mode {}",
                    static_cast<int>(mode)));
        }
        const auto& X = bayes::get_matrix_ref(it->X);
        stats.push_back({
            .marker_variance_sum = stats::detail::var(X).sum(),
            .heritability_init
            = heritability_for(mode) * model.phenotype_variance(),
        });
    }
    return stats;
}

}  // namespace gelex

namespace gelex::bayes
{

auto build_bayes_method(
    const BayesConfig& config,
    std::span<const GeneticStats> stats,
    double phenotype_variance) -> BayesMethod
{
    const auto& policy = policy_for(config.base);

    GeneticPrior prior;
    switch (config.mode)
    {
        case GeneticMode::A:
            if (stats.size() != 1)
            {
                throw GelexException(
                    "Single-effect mode requires one GeneticStats entry");
            }
            prior.spec = GeneticSpec::make(stats[0], policy);
            break;
        case GeneticMode::D:
            if (stats.size() != 1)
            {
                throw GelexException(
                    "Single-effect mode requires one GeneticStats entry");
            }
            prior.spec = GeneticSpec::make(stats[0], policy, config.dominance);
            break;
        case GeneticMode::AD:
            if (stats.size() != 2)
            {
                throw GelexException(
                    "A+D mode requires two GeneticStats entries");
            }
            prior.spec = JointSpec{
                GeneticSpec::make(stats[0], policy),
                GeneticSpec::make(stats[1], policy, config.dominance),
            };
            break;
        case GeneticMode::kCount:
            throw GelexException("Invalid GeneticMode: kCount");
    }
    prior.mixture = Mixture::make(policy, config.estimate_pi);

    BayesMethod method;
    method.config = config;
    method.genetics.push_back(std::move(prior));
    method.residual = VarianceSpec::make(phenotype_variance);
    return method;
}

}  // namespace gelex::bayes
