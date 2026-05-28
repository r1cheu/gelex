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
#include <random>
#include <vector>

#include <Eigen/Core>
#include <catch2/catch_test_macros.hpp>

#include "gelex/exception.h"
#include "gelex/io/mcmc/checkpoint_reader.h"
#include "gelex/io/mcmc/checkpoint_writer.h"
#include "gelex/model/bayes/genetic_priors/gaussian.h"
#include "gelex/model/bayes/model.h"
#include "gelex/model/bayes/prior.h"
#include "gelex/model/bayes/prior_parameters.h"
#include "gelex/model/bayes/state.h"
#include "gelex/types/fixed_designs.h"
#include "gelex/types/genetic_effect_type.h"
#include "genotype_fixture.h"

namespace
{

auto make_variance(double value) -> gelex::bayes::VarianceParameter
{
    return gelex::bayes::VarianceParameter{
        value, gelex::bayes::ScaledInvChiSqPrior{4.0, 1.0}};
}

auto make_genotype() -> gelex::genotype::Genotype
{
    Eigen::MatrixXd data{{0.0, 1.0}, {1.0, 0.0}, {2.0, 1.0}};
    auto mean = data.colwise().mean().transpose().eval();
    auto stddev = Eigen::VectorXd::Ones(data.cols());
    return gelex::test::GenotypeBuilder::build(
        std::move(data), std::move(mean), std::move(stddev));
}

auto make_model() -> gelex::BayesModel
{
    std::vector<gelex::bayes::GeneticDesign> genetics;
    genetics.emplace_back(gelex::GeneticMode::A, make_genotype());

    return gelex::BayesModel{
        Eigen::VectorXd{{1.0, 2.0, 3.0}},
        gelex::FixedDesign::build(3),
        {},
        std::move(genetics)};
}

auto make_prior() -> gelex::bayes::BayesPrior
{
    std::vector<gelex::bayes::GeneticPriorBlock> genetics;
    genetics.emplace_back(
        std::make_unique<gelex::bayes::SingleSharedGaussianPrior>(
            gelex::GeneticMode::A,
            gelex::bayes::SharedMarkerVariance{make_variance(0.1)}));

    return gelex::bayes::BayesPrior{
        gelex::bayes::RandomPrior{make_variance(0.3)},
        std::move(genetics),
        gelex::bayes::ResidualPrior{make_variance(0.4)}};
}

auto make_state() -> gelex::BayesState
{
    auto model = make_model();
    auto prior = make_prior();
    return gelex::BayesState{model, prior};
}

}  // namespace

TEST_CASE("MCMC checkpoint writer is not implemented", "[checkpoint]")
{
    auto state = make_state();
    std::mt19937_64 rng{123};

    REQUIRE_THROWS_AS(
        gelex::write_checkpoint(state, rng, "/tmp/gelex_unused"),
        gelex::GelexException);
}

TEST_CASE("MCMC checkpoint reader is not implemented", "[checkpoint]")
{
    auto state = make_state();

    REQUIRE_THROWS_AS(
        gelex::read_checkpoint("/tmp/gelex_unused.ckpt", state),
        gelex::GelexException);
}
