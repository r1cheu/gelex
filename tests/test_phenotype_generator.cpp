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
#include <random>

#include <Eigen/Core>
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include "gelex/data/genetic_value_calculator.h"
#include "gelex/data/phenotype_generator.h"
#include "utils/math_utils.h"

using namespace gelex;  // NOLINT
using Catch::Matchers::WithinAbs;

namespace
{
constexpr double VARIANCE_TOLERANCE = 0.15;
}

TEST_CASE("PhenotypeGenerator - additive only", "[phenotype_generator]")
{
    constexpr Eigen::Index N_SAMPLES = 1000;
    constexpr double H2 = 0.5;

    std::mt19937_64 setup_rng(42);
    Eigen::VectorXd additive_values(N_SAMPLES);
    for (Eigen::Index i = 0; i < N_SAMPLES; ++i)
    {
        additive_values(i) = std::normal_distribution<double>(0.0, 1.0)(setup_rng);
    }

    GeneticValues gv{
        .additive = additive_values,
        .dominance = Eigen::VectorXd::Zero(N_SAMPLES),
    };

    PhenotypeGeneratorConfig config{
        .h2 = H2, .d2 = 0.0, .intercept = 0.0, .seed = 123};
    PhenotypeGenerator generator(config);

    auto result = generator.generate(gv);

    SECTION("Output size matches input")
    {
        REQUIRE(result.phenotypes.size() == N_SAMPLES);
    }

    SECTION("True h2 is close to target")
    {
        REQUIRE_THAT(result.true_h2, WithinAbs(H2, VARIANCE_TOLERANCE));
    }

    SECTION("True d2 is zero when d2 config is zero")
    {
        REQUIRE(result.true_d2 == 0.0);
    }

    SECTION("Dominance unchanged when no dominance")
    {
        REQUIRE(gv.dominance.isZero());
    }
}

TEST_CASE("PhenotypeGenerator - with dominance", "[phenotype_generator]")
{
    constexpr Eigen::Index N_SAMPLES = 1000;
    constexpr double H2 = 0.4;
    constexpr double D2 = 0.2;

    std::mt19937_64 setup_rng(42);
    Eigen::VectorXd additive_values(N_SAMPLES);
    Eigen::VectorXd dominance_values(N_SAMPLES);
    for (Eigen::Index i = 0; i < N_SAMPLES; ++i)
    {
        additive_values(i) = std::normal_distribution<double>(0.0, 1.0)(setup_rng);
        dominance_values(i) = std::normal_distribution<double>(0.0, 0.5)(setup_rng);
    }

    GeneticValues gv{
        .additive = additive_values,
        .dominance = dominance_values,
    };

    PhenotypeGeneratorConfig config{
        .h2 = H2, .d2 = D2, .intercept = 0.0, .seed = 123};
    PhenotypeGenerator generator(config);

    auto result = generator.generate(gv);

    SECTION("True h2 is close to target")
    {
        REQUIRE_THAT(result.true_h2, WithinAbs(H2, VARIANCE_TOLERANCE));
    }

    SECTION("True d2 is close to target")
    {
        REQUIRE_THAT(result.true_d2, WithinAbs(D2, VARIANCE_TOLERANCE));
    }

    SECTION("Dominance values are scaled in place")
    {
        double genetic_var = detail::var(gv.additive)(0);
        double dom_var = detail::var(gv.dominance)(0);
        double target_dom_var = genetic_var * D2 / H2;
        REQUIRE_THAT(dom_var, WithinAbs(target_dom_var, target_dom_var * VARIANCE_TOLERANCE));
    }
}

TEST_CASE("PhenotypeGenerator - intercept", "[phenotype_generator]")
{
    constexpr Eigen::Index N_SAMPLES = 100;
    constexpr double INTERCEPT = 10.0;

    Eigen::VectorXd additive_values = Eigen::VectorXd::Zero(N_SAMPLES);

    GeneticValues gv{
        .additive = additive_values,
        .dominance = Eigen::VectorXd::Zero(N_SAMPLES),
    };

    PhenotypeGeneratorConfig config{
        .h2 = 0.5, .d2 = 0.0, .intercept = INTERCEPT, .seed = 123};
    PhenotypeGenerator generator(config);

    auto result = generator.generate(gv);

    double mean_phenotype = result.phenotypes.mean();
    REQUIRE_THAT(mean_phenotype, WithinAbs(INTERCEPT, 1.0));
}

TEST_CASE("PhenotypeGenerator - reproducibility", "[phenotype_generator]")
{
    constexpr Eigen::Index N_SAMPLES = 100;

    std::mt19937_64 setup_rng(42);
    Eigen::VectorXd additive_values(N_SAMPLES);
    for (Eigen::Index i = 0; i < N_SAMPLES; ++i)
    {
        additive_values(i) = std::normal_distribution<double>(0.0, 1.0)(setup_rng);
    }

    GeneticValues gv{
        .additive = additive_values,
        .dominance = Eigen::VectorXd::Zero(N_SAMPLES),
    };

    PhenotypeGeneratorConfig config{
        .h2 = 0.5, .d2 = 0.0, .intercept = 0.0, .seed = 123};
    PhenotypeGenerator generator1(config);
    auto result1 = generator1.generate(gv);

    PhenotypeGenerator generator2(config);
    auto result2 = generator2.generate(gv);

    for (Eigen::Index i = 0; i < N_SAMPLES; ++i)
    {
        REQUIRE(result1.phenotypes(i) == result2.phenotypes(i));
    }
}
