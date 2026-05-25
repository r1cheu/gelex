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

#include <cmath>
#include <string>
#include <utility>
#include <vector>

#include <Eigen/Core>
#include <catch2/catch_test_macros.hpp>

#include "gelex/data/genotype/genotype.h"
#include "gelex/exception.h"
#include "gelex/infra/stats/detail/var.h"
#include "gelex/model/bayes/effects.h"
#include "gelex/model/bayes/model.h"
#include "gelex/types/fixed_effects.h"
#include "gelex/types/genetic_effect_type.h"
#include "genotype_fixture.h"

using gelex::BayesModel;
using gelex::FixedEffect;
using gelex::GelexException;
using gelex::GeneticMode;

namespace
{

auto make_genotype(Eigen::MatrixXd data) -> gelex::genotype::Genotype
{
    const Eigen::Index cols = data.cols();
    auto mean = data.colwise().mean().transpose().eval();
    auto stddev = Eigen::VectorXd::Ones(cols);
    return gelex::test::GenotypeBuilder::build(
        std::move(data), std::move(mean), std::move(stddev));
}

auto make_genetic_effect(GeneticMode mode, Eigen::MatrixXd data)
    -> gelex::bayes::GeneticEffect
{
    return gelex::bayes::GeneticEffect{mode, make_genotype(std::move(data))};
}

auto make_phenotype() -> Eigen::VectorXd
{
    return Eigen::VectorXd{{1.0, 2.0, 4.0}};
}

auto make_genetic_data(Eigen::Index rows = 3) -> Eigen::MatrixXd
{
    Eigen::MatrixXd data(rows, 2);
    for (Eigen::Index i = 0; i < rows; ++i)
    {
        data(i, 0) = static_cast<double>(i);
        data(i, 1) = static_cast<double>(i + 1);
    }
    return data;
}

}  // namespace

TEST_CASE("BayesModel rejects effect row mismatches", "[bayes_model]")
{
    SECTION("fixed")
    {
        std::vector<gelex::bayes::GeneticEffect> genetics;
        genetics.push_back(
            make_genetic_effect(GeneticMode::A, make_genetic_data()));

        REQUIRE_THROWS_AS(
            BayesModel(
                make_phenotype(),
                FixedEffect::build(2),
                {},
                std::move(genetics)),
            GelexException);
    }

    SECTION("random")
    {
        std::vector<gelex::bayes::RandomEffect> random;
        random.emplace_back(
            "batch",
            std::vector<std::string>{"a", "b"},
            Eigen::MatrixXd{{1.0, 0.0}, {0.0, 1.0}});

        std::vector<gelex::bayes::GeneticEffect> genetics;
        genetics.push_back(
            make_genetic_effect(GeneticMode::A, make_genetic_data()));

        REQUIRE_THROWS_AS(
            BayesModel(
                make_phenotype(),
                FixedEffect::build(3),
                std::move(random),
                std::move(genetics)),
            GelexException);
    }

    SECTION("genetic")
    {
        std::vector<gelex::bayes::GeneticEffect> genetics;
        genetics.push_back(
            make_genetic_effect(GeneticMode::A, make_genetic_data(2)));

        REQUIRE_THROWS_AS(
            BayesModel(
                make_phenotype(),
                FixedEffect::build(3),
                {},
                std::move(genetics)),
            GelexException);
    }
}

TEST_CASE("BayesModel rejects duplicate genetic modes", "[bayes_model]")
{
    std::vector<gelex::bayes::GeneticEffect> genetics;
    genetics.push_back(
        make_genetic_effect(GeneticMode::A, make_genetic_data()));
    genetics.push_back(
        make_genetic_effect(GeneticMode::A, make_genetic_data()));

    REQUIRE_THROWS_AS(
        BayesModel(
            make_phenotype(), FixedEffect::build(3), {}, std::move(genetics)),
        GelexException);
}

TEST_CASE(
    "GeneticEffect caches design_variance from design matrix",
    "[bayes_model]")
{
    const Eigen::MatrixXd data{{0.0, 1.0}, {1.0, 1.5}, {2.0, 3.0}};
    const double expected = gelex::stats::detail::var<0>(data).sum();

    const auto effect = make_genetic_effect(GeneticMode::A, data);

    REQUIRE(std::abs(effect.design_variance - expected) < 1e-12);
}
