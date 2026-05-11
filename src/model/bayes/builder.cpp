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
#include <variant>
#include <vector>

#include "gelex/data/genotype/storage.h"
#include "gelex/data/pipe/geno.h"
#include "gelex/data/pipe/pheno.h"
#include "gelex/infra/stats/descriptive.h"
#include "gelex/model/bayes/algorithm_shape.h"
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
    const bayes::BayesConfig& /*config*/) -> std::vector<bayes::GeneticStats>
{
    constexpr double kAdditiveHeritability = 0.5;
    constexpr double kDominanceHeritability = 0.2;

    std::vector<bayes::GeneticStats> stats;
    stats.reserve(model.genetics().size());
    for (const auto& effect : model.genetics())
    {
        const double h2 = (effect.type == GeneticMode::D)
                              ? kDominanceHeritability
                              : kAdditiveHeritability;
        const auto& X = bayes::get_matrix_ref(effect.X);
        stats.push_back({
            .mode = effect.type,
            .marker_variance_sum = stats::detail::var(X).sum(),
            .heritability_init = h2 * model.phenotype_variance(),
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

    std::vector<GeneticMode> requested;
    requested.reserve(stats.size());
    for (const auto& s : stats)
    {
        requested.push_back(s.mode);
    }
    const auto shape = resolve_shape(policy, requested);

    auto make_prior = [&](std::variant<GeneticSpec, JointSpec> spec)
    {
        GeneticPrior prior;
        prior.spec = std::move(spec);
        prior.mixture = Mixture::make(policy, config.estimate_pi);
        return prior;
    };

    auto make_spec = [&](const GeneticStats& s) -> GeneticSpec
    {
        return s.mode == GeneticMode::D
                   ? GeneticSpec::make(s, policy, config.dominance)
                   : GeneticSpec::make(s, policy);
    };

    BayesMethod method;
    method.config = config;
    method.residual = VarianceSpec::make(phenotype_variance);

    if (shape == AlgorithmShape::ad_joint)
    {
        method.genetics.push_back(
            make_prior(JointSpec{make_spec(stats[0]), make_spec(stats[1])}));
    }
    else
    {
        for (const auto& s : stats)
        {
            method.genetics.push_back(make_prior(make_spec(s)));
        }
    }
    return method;
}

}  // namespace gelex::bayes
