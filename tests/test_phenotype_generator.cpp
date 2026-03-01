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

#include "gelex/algo/sim/genetic_value_calculator.h"
#include "gelex/algo/sim/phenotype_generator.h"
#include "gelex/infra/utils/math_utils.h"

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

    std::mt19937_64 rng(123);
    PhenotypeGenerator generator(H2, 0.0, 0.0, rng);
    std::optional<double> observed_h2;
    std::optional<double> observed_d2;
    auto phenotypes = generator.generate(gv, [&](const SimulateEvent& event) {
        if (const auto* h2_event = std::get_if<HeritabilityGeneratedEvent>(&event))
        {
            observed_h2 = h2_event->additive;
            observed_d2 = h2_event->dominance;
        }
    });

    SECTION("Output size matches input")
    {
        REQUIRE(phenotypes.size() == N_SAMPLES);
    }

    SECTION("True h2 is close to target")
    {
        REQUIRE(observed_h2.has_value());
        REQUIRE_THAT(*observed_h2, WithinAbs(H2, VARIANCE_TOLERANCE));
    }

    SECTION("True d2 is zero when d2 config is zero")
    {
        REQUIRE(observed_d2 == std::nullopt);
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

    std::mt19937_64 rng(123);
    PhenotypeGenerator generator(H2, D2, 0.0, rng);
    std::optional<double> observed_h2;
    std::optional<double> observed_d2;
    auto phenotypes = generator.generate(gv, [&](const SimulateEvent& event) {
        if (const auto* h2_event = std::get_if<HeritabilityGeneratedEvent>(&event))
        {
            observed_h2 = h2_event->additive;
            observed_d2 = h2_event->dominance;
        }
    });

    SECTION("True h2 is close to target")
    {
        (void)phenotypes;
        REQUIRE(observed_h2.has_value());
        REQUIRE_THAT(*observed_h2, WithinAbs(H2, VARIANCE_TOLERANCE));
    }

    SECTION("True d2 is close to target")
    {
        REQUIRE(observed_d2.has_value());
        REQUIRE_THAT(*observed_d2, WithinAbs(D2, VARIANCE_TOLERANCE));
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

    std::mt19937_64 rng(123);
    PhenotypeGenerator generator(0.5, 0.0, INTERCEPT, rng);

    auto phenotypes = generator.generate(gv);

    double mean_phenotype = phenotypes.mean();
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

    std::mt19937_64 rng1(123);
    PhenotypeGenerator generator1(0.5, 0.0, 0.0, rng1);
    auto result1 = generator1.generate(gv);

    std::mt19937_64 rng2(123);
    PhenotypeGenerator generator2(0.5, 0.0, 0.0, rng2);
    auto result2 = generator2.generate(gv);

    for (Eigen::Index i = 0; i < N_SAMPLES; ++i)
    {
        REQUIRE(result1(i) == result2(i));
    }
}
