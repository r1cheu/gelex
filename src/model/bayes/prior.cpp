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

#include "gelex/model/bayes/prior.h"

#include <algorithm>
#include <cstddef>
#include <optional>
#include <span>
#include <string_view>
#include <utility>

#include <fmt/format.h>
#include <Eigen/Core>

#include "gelex/exception.h"
#include "gelex/model/bayes/bayes_policy.h"
#include "gelex/model/bayes/constrain_vector.h"
#include "gelex/model/bayes/prior_constants.h"

namespace gelex::bayes
{

MarkerVarianceSpec::MarkerVarianceSpec(
    MarkerVarianceScope scope,
    VarianceSpec variance)
    : scope_(scope), variance_(variance)
{
}

ProportionSpec::ProportionSpec(
    Simplex<double> initial_value,
    DirichletPrior prior,
    ProportionUpdate update)
    : initial_value_(std::move(initial_value)),
      prior_(std::move(prior)),
      update_(update)
{
    if (initial_value_.size() != prior_.concentration.size())
    {
        throw GelexException(
            "ProportionSpec: initial value and prior concentration must have "
            "the same size");
    }
}

auto OldVarianceSpec::make(double phenotype_variance) -> OldVarianceSpec
{
    constexpr double kResidualVarianceProportion = 0.3;
    return OldVarianceSpec{
        .scope = MarkerVarianceScope::per_effect,
        .init = kResidualVarianceProportion * phenotype_variance,
        .prior = ScaledInvChiSqPrior{
            prior_constants::RESIDUAL_VARIANCE_SHAPE,
            prior_constants::RESIDUAL_VARIANCE_SCALE,
        },
    };
}

auto OldVarianceSpec::make(double init, MarkerVarianceScope scope)
    -> OldVarianceSpec
{
    return OldVarianceSpec{
        .scope = scope,
        .init = init,
        .prior = ScaledInvChiSqPrior{
            prior_constants::MARKER_VARIANCE_SHAPE,
            prior_constants::MARKER_VARIANCE_SCALE_MULTIPLIER * init,
        },
    };
}

auto Mixture::make(const BayesPolicy& policy, bool estimate_pi)
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
            OldDirichletPrior{
                Eigen::VectorXi::Ones(m.default_proportion.size())},
            estimate_pi,
        },
    };
}

BayesPrior::BayesPrior(
    VarianceSpec random,
    std::vector<std::unique_ptr<GeneticPrior>> genetics,
    VarianceSpec residual)
    : random_(random), genetics_(std::move(genetics)), residual_(residual)
{
    const auto check_variance
        = [](const VarianceSpec& spec, std::string_view context)
    {
        if (spec.initial_value <= 0.0)
        {
            throw GelexException(
                fmt::format(
                    "BayesPrior: {} variance initial_value must be > 0",
                    context));
        }
        if (spec.prior.degrees_of_freedom <= 0.0)
        {
            throw GelexException(
                fmt::format(
                    "BayesPrior: {} variance prior degrees_of_freedom must be "
                    "> 0",
                    context));
        }
        if (spec.prior.scale <= 0.0)
        {
            throw GelexException(
                fmt::format(
                    "BayesPrior: {} variance prior scale must be > 0",
                    context));
        }
    };

    check_variance(random_, "random");
    check_variance(residual_, "residual");

    std::vector<GeneticMode> seen_modes;
    for (const auto& block : genetics_)
    {
        if (block == nullptr)
        {
            throw GelexException("BayesPrior: null GeneticPrior block");
        }
        const auto block_modes = block->modes();
        if (block_modes.empty())
        {
            throw GelexException("BayesPrior: GeneticPrior block has no modes");
        }
        for (const auto mode : block_modes)
        {
            if (std::ranges::contains(seen_modes, mode))
            {
                throw GelexException(
                    fmt::format(
                        "BayesPrior: duplicate GeneticMode {} across blocks",
                        mode));
            }
            seen_modes.push_back(mode);
        }
    }
}

}  // namespace gelex::bayes
