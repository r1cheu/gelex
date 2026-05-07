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
#include <optional>
#include <utility>
#include <variant>
#include <vector>

#include <fmt/format.h>

#include <Eigen/Core>

#include "gelex/data/genotype/storage.h"
#include "gelex/data/pipe/geno.h"
#include "gelex/data/pipe/pheno.h"
#include "gelex/exception.h"
#include "gelex/infra/stats/conjugate_prior.h"
#include "gelex/infra/stats/descriptive.h"
#include "gelex/model/bayes/bayes_policy.h"
#include "gelex/model/bayes/effects.h"
#include "gelex/model/bayes/method.h"
#include "gelex/model/bayes/model.h"
#include "gelex/model/bayes/prior.h"
#include "gelex/model/bayes/prior_constants.h"
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

namespace
{

constexpr double kAdditiveHeritability = 0.5;
constexpr double kDominanceHeritability = 0.2;
constexpr double kResidualVarianceProportion = 0.3;

auto heritability_for(GeneticMode mode) -> double
{
    return mode == GeneticMode::D ? kDominanceHeritability
                                  : kAdditiveHeritability;
}

auto compute_genetic_stats(
    const std::vector<bayes::GeneticEffect>& genetics,
    GeneticMode mode,
    double phenotype_variance) -> bayes::GeneticStats
{
    auto it = std::ranges::find(genetics, mode, &bayes::GeneticEffect::type);
    if (it == genetics.end())
    {
        throw GelexException(
            fmt::format(
                "BayesModel missing genetic effect for mode {}",
                static_cast<int>(mode)));
    }
    const auto& X = bayes::get_matrix_ref(it->X);
    return {
        .marker_variance_sum = stats::detail::var(X).sum(),
        .heritability_init = heritability_for(mode) * phenotype_variance,
    };
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

void apply_proportion_override(
    bayes::BayesMethod& method,
    GeneticMode mode,
    const std::optional<std::vector<double>>& values)
{
    if (!values)
    {
        return;
    }
    for (auto& prior : method.genetics)
    {
        if (find_spec(prior, mode) != nullptr && prior.mixture)
        {
            prior.mixture->proportions.init = Eigen::Map<const Eigen::VectorXd>(
                values->data(), static_cast<Eigen::Index>(values->size()));
        }
    }
}

void apply_multiplier_override(
    bayes::BayesMethod& method,
    GeneticMode mode,
    const std::optional<std::vector<double>>& values)
{
    if (!values)
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
            sm->multiplier = Eigen::Map<const Eigen::VectorXd>(
                values->data(), static_cast<Eigen::Index>(values->size()));
        }
    }
}

void apply_positive_prob_override(bayes::BayesMethod& method, double value)
{
    for (auto& prior : method.genetics)
    {
        if (auto* spec = find_spec(prior, GeneticMode::D);
            spec != nullptr && spec->sign)
        {
            spec->sign->init = Eigen::VectorXd{{value, 1.0 - value}};
        }
    }
}

auto residual_spec(double phenotype_variance) -> bayes::VarianceSpec
{
    return bayes::VarianceSpec{
        bayes::VarianceScope::per_block,
        kResidualVarianceProportion * phenotype_variance,
        bayes::ScaledInvChiSqPrior{
            prior_constants::RESIDUAL_VARIANCE_SHAPE,
            prior_constants::RESIDUAL_VARIANCE_SCALE,
        },
    };
}

}  // namespace

auto build_bayes_method(
    const PriorOverrides& overrides,
    const BayesModel& model) -> bayes::BayesMethod
{
    const auto& genetics = model.genetics();
    const double phenotype_variance = overrides.phenotype_variance;

    std::vector<bayes::GeneticStats> stats;
    switch (overrides.method.mode)
    {
        case GeneticMode::A:
            stats.push_back(compute_genetic_stats(
                genetics, GeneticMode::A, phenotype_variance));
            break;
        case GeneticMode::D:
            stats.push_back(compute_genetic_stats(
                genetics, GeneticMode::D, phenotype_variance));
            break;
        case GeneticMode::AD:
            stats.push_back(compute_genetic_stats(
                genetics, GeneticMode::A, phenotype_variance));
            stats.push_back(compute_genetic_stats(
                genetics, GeneticMode::D, phenotype_variance));
            break;
    }

    auto method = bayes::build_bayes_method(overrides.method, stats);
    apply_proportion_override(method, GeneticMode::A, overrides.pi);
    apply_proportion_override(method, GeneticMode::D, overrides.dpi);
    apply_multiplier_override(method, GeneticMode::A, overrides.multiplier);
    apply_multiplier_override(method, GeneticMode::D, overrides.dmultiplier);
    apply_positive_prob_override(method, overrides.positive_prob);
    method.residual = residual_spec(phenotype_variance);

    return method;
}

}  // namespace gelex

namespace gelex::bayes
{

namespace
{

auto make_variance_prior(double init) -> ScaledInvChiSqPrior
{
    return {
        prior_constants::MARKER_VARIANCE_SHAPE,
        prior_constants::MARKER_VARIANCE_SCALE_MULTIPLIER * init,
    };
}

auto build_spec(
    const GeneticStats& stats,
    const BayesPolicy& policy,
    double non_zero,
    GeneticMode mode,
    DominancePolicy dominance) -> GeneticSpec
{
    if (stats.heritability_init <= 0.0)
    {
        throw GelexException("Block target variance must be positive");
    }
    const double num_non_zero = stats.marker_variance_sum * non_zero;
    if (num_non_zero <= 0.0)
    {
        throw GelexException("Number of non-zero SNPs must be positive");
    }
    const double init = stats.heritability_init / num_non_zero;

    GeneticSpec spec;
    spec.mode = mode;
    spec.variance
        = VarianceSpec{policy.variance_scope, init, make_variance_prior(init)};

    if (dominance == DominancePolicy::asymmetric && mode == GeneticMode::D)
    {
        spec.sign = CategoricalSpec{
            Eigen::VectorXd{{0.5, 0.5}},
            DirichletPrior{Eigen::VectorXi{{1, 1}}},
            false,
        };
    }
    return spec;
}

auto build_mixture(const BayesPolicy& policy, bool estimate_pi)
    -> std::optional<Mixture>
{
    if (!policy.mixture)
    {
        return std::nullopt;
    }
    const auto& m = *policy.mixture;
    return Mixture{
        m.strategy,
        CategoricalSpec{
            m.default_proportion,
            DirichletPrior{Eigen::VectorXi::Ones(m.default_proportion.size())},
            estimate_pi,
        },
    };
}

}  // namespace

auto build_bayes_method(
    const BayesConfig& config,
    std::span<const GeneticStats> stats) -> BayesMethod
{
    const auto& policy = policy_for(config.base);
    const double non_zero = policy.mixture
                                ? 1.0 - policy.mixture->default_proportion[0]
                                : prior_constants::NON_MIXTURE_PROPORTION;

    GeneticPrior prior;
    if (config.mode == GeneticMode::AD)
    {
        if (stats.size() != 2)
        {
            throw GelexException("A+D mode requires two GeneticStats entries");
        }
        prior.spec = JointSpec{
            build_spec(
                stats[0], policy, non_zero, GeneticMode::A, config.dominance),
            build_spec(
                stats[1], policy, non_zero, GeneticMode::D, config.dominance),
        };
    }
    else
    {
        if (stats.size() != 1)
        {
            throw GelexException(
                "Single-effect mode requires one GeneticStats entry");
        }
        prior.spec = build_spec(
            stats[0], policy, non_zero, config.mode, config.dominance);
    }
    prior.mixture = build_mixture(policy, config.estimate_pi);

    BayesMethod method;
    method.config = config;
    method.genetics.push_back(std::move(prior));
    return method;
}

}  // namespace gelex::bayes
