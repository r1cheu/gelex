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

#include "gelex/model/bayes/builder_.h"

#include <utility>

#include <Eigen/Core>

#include "gelex/exception.h"
#include "gelex/model/bayes/bayes_policy.h"
#include "gelex/model/bayes/method_.h"
#include "gelex/model/bayes/prior_.h"
#include "gelex/model/bayes/prior_constants.h"
#include "gelex/types/genetic_effect_type.h"

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

    GeneticTerm term;
    if (config.mode == GeneticMode::AD)
    {
        if (stats.size() != 2)
        {
            throw GelexException("A+D mode requires two GeneticStats entries");
        }
        term.spec = JointSpec{
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
        term.spec = build_spec(
            stats[0], policy, non_zero, config.mode, config.dominance);
    }
    term.mixture = build_mixture(policy, config.estimate_pi);

    BayesMethod method;
    method.config = config;
    method.genetic_terms.push_back(std::move(term));
    return method;
}

}  // namespace gelex::bayes
