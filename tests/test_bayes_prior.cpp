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

#include <memory>
#include <utility>
#include <variant>
#include <vector>

#include <catch2/catch_test_macros.hpp>

#include "gelex/exception.h"
#include "gelex/model/bayes/capabilities.h"
#include "gelex/model/bayes/genetic_prior.h"
#include "gelex/model/bayes/genetic_priors/gaussian.h"
#include "gelex/model/bayes/prior.h"
#include "gelex/model/bayes/prior_parameters.h"
#include "gelex/model/bayes/state_capabilities.h"
#include "gelex/types/genetic_effect_type.h"

namespace
{

auto make_variance(double value) -> gelex::bayes::VarianceParameter
{
    return gelex::bayes::VarianceParameter{
        value, gelex::bayes::ScaledInvChiSqPrior{4.0, 1.0}};
}

auto make_random_prior() -> gelex::bayes::RandomPrior
{
    return gelex::bayes::RandomPrior{make_variance(0.2)};
}

auto make_residual_prior() -> gelex::bayes::ResidualPrior
{
    return gelex::bayes::ResidualPrior{make_variance(0.8)};
}

auto make_shared_prior(gelex::GeneticMode mode)
    -> std::unique_ptr<gelex::bayes::SingleGeneticPrior>
{
    return std::make_unique<gelex::bayes::SingleSharedGaussianPrior>(
        mode, gelex::bayes::SharedMarkerVariance{make_variance(0.1)});
}

}  // namespace

TEST_CASE("BayesPrior owns single and joint genetic blocks", "[bayes_prior]")
{
    std::vector<gelex::bayes::GeneticPriorBlock> genetics;
    genetics.emplace_back(make_shared_prior(gelex::GeneticMode::A));

    gelex::bayes::BayesPrior prior{
        make_random_prior(), std::move(genetics), make_residual_prior()};

    REQUIRE(prior.genetics().size() == 1);
    const auto* single
        = std::get_if<std::unique_ptr<gelex::bayes::SingleGeneticPrior>>(
            &prior.genetics()[0]);
    REQUIRE(single != nullptr);
    REQUIRE((*single)->mode() == gelex::GeneticMode::A);
}

TEST_CASE("Genetic prior get and get_if expose capabilities", "[bayes_prior]")
{
    gelex::bayes::SingleSharedGaussianPrior prior{
        gelex::GeneticMode::A,
        gelex::bayes::SharedMarkerVariance{make_variance(0.1)}};

    REQUIRE(prior.get_if<gelex::bayes::SingleSharedMarkerVarCap>() != nullptr);
    REQUIRE(
        prior.get_if<gelex::bayes::SingleFixedMixtureProportionCap>()
        == nullptr);
    REQUIRE_THROWS_AS(
        prior.get<gelex::bayes::SingleFixedMixtureProportionCap>(),
        gelex::GelexException);
}

TEST_CASE("Genetic prior states expose capabilities", "[bayes_prior]")
{
    gelex::bayes::SingleSharedGaussianPrior prior{
        gelex::GeneticMode::A,
        gelex::bayes::SharedMarkerVariance{make_variance(0.1)}};
    auto state = prior.make_state(3, 2);

    auto& variance
        = state->get<gelex::bayes::SingleSharedVarianceStateCap>().variance();
    REQUIRE(variance == 0.1);
    REQUIRE(
        state->get_if<gelex::bayes::SingleMixtureAssignmentStateCap>()
        == nullptr);
}

TEST_CASE("BayesPrior rejects duplicate genetic modes", "[bayes_prior]")
{
    std::vector<gelex::bayes::GeneticPriorBlock> genetics;
    genetics.emplace_back(make_shared_prior(gelex::GeneticMode::A));
    genetics.emplace_back(make_shared_prior(gelex::GeneticMode::A));

    REQUIRE_THROWS_AS(
        gelex::bayes::BayesPrior(
            make_random_prior(), std::move(genetics), make_residual_prior()),
        gelex::GelexException);
}
