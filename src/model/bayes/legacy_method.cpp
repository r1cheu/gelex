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

#include "gelex/model/bayes/legacy_method.h"

#include <optional>

#include <Eigen/Core>

#include "gelex/exception.h"
#include "gelex/model/bayes/legacy_bayes_policy.h"
#include "gelex/model/bayes/legacy_builder.h"
#include "gelex/model/bayes/legacy_prior_constants.h"
#include "gelex/model/bayes/prior.h"

namespace gelex::bayes
{

auto GeneticSpec::make(const GeneticStats& stats, const BayesPolicy& policy)
    -> GeneticSpec
{
    if (stats.heritability_init <= 0.0)
    {
        throw GelexException("Block target variance must be positive");
    }
    const double non_zero = policy.mixture
                                ? 1.0 - policy.mixture->default_proportion[0]
                                : prior_constants::NON_MIXTURE_PROPORTION;
    const double num_non_zero = stats.marker_variance_sum * non_zero;
    if (num_non_zero <= 0.0)
    {
        throw GelexException("Number of non-zero SNPs must be positive");
    }
    const double init = stats.heritability_init / num_non_zero;

    return GeneticSpec{
        .mode = GeneticMode::A,
        .variance = OldVarianceSpec::make(init, policy.variance_scope),
        .sign = std::nullopt,
    };
}

auto GeneticSpec::make(
    const GeneticStats& stats,
    const BayesPolicy& policy,
    DominancePolicy dominance) -> GeneticSpec
{
    auto spec = make(stats, policy);
    spec.mode = GeneticMode::D;
    if (dominance == DominancePolicy::asymmetric)
    {
        spec.sign = CategoricalSpec{
            Eigen::VectorXd{{0.5, 0.5}},
            OldDirichletPrior{Eigen::VectorXi{{1, 1}}},
            false,
        };
    }
    return spec;
}

}  // namespace gelex::bayes
