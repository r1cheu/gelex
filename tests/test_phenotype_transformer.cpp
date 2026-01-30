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

#include <Eigen/Dense>
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include "utils/phenotype_transformer.h"

using gelex::detail::PhenotypeTransformer;

TEST_CASE(
    "PhenotypeTransformer - DINT produces approximately standard normal",
    "[phenotype_transformer]")
{
    std::mt19937 gen(12345);
    std::lognormal_distribution<double> dist(0.0, 1.0);

    Eigen::VectorXd phenotype(1000);
    for (int i = 0; i < 1000; ++i)
    {
        phenotype[i] = dist(gen);
    }

    PhenotypeTransformer transformer;
    transformer.apply_dint(phenotype);

    auto mean = phenotype.mean();
    auto variance = (phenotype.array() - mean).square().mean();
    auto std_dev = std::sqrt(variance);

    REQUIRE_THAT(mean, Catch::Matchers::WithinAbs(0.0, 0.1));
    REQUIRE_THAT(std_dev, Catch::Matchers::WithinAbs(1.0, 0.1));
}

TEST_CASE(
    "PhenotypeTransformer - DINT handles ties correctly",
    "[phenotype_transformer]")
{
    Eigen::VectorXd phenotype{{1.0, 2.0, 2.0, 2.0, 3.0, 3.0, 4.0, 5.0}};

    PhenotypeTransformer transformer;
    transformer.apply_dint(phenotype);

    REQUIRE(std::abs(phenotype[1] - phenotype[2]) < 1e-10);
    REQUIRE(std::abs(phenotype[2] - phenotype[3]) < 1e-10);
    REQUIRE(std::abs(phenotype[4] - phenotype[5]) < 1e-10);

    REQUIRE(phenotype[0] < phenotype[1]);
    REQUIRE(phenotype[3] < phenotype[4]);
    REQUIRE(phenotype[5] < phenotype[6]);
    REQUIRE(phenotype[6] < phenotype[7]);
}

TEST_CASE(
    "PhenotypeTransformer - IINT with covariates",
    "[phenotype_transformer]")
{
    std::mt19937 gen(42);
    std::normal_distribution<double> normal_dist(0.0, 1.0);
    std::lognormal_distribution<double> lognormal_dist(0.0, 0.5);

    int n = 500;

    Eigen::MatrixXd covariates(n, 3);
    covariates.col(0).setOnes();
    for (int i = 0; i < n; ++i)
    {
        covariates(i, 1) = normal_dist(gen);
        covariates(i, 2) = normal_dist(gen);
    }

    Eigen::VectorXd beta{{2.0, 0.5, -0.3}};

    Eigen::VectorXd phenotype(n);
    for (int i = 0; i < n; ++i)
    {
        phenotype[i] = covariates.row(i).dot(beta) + lognormal_dist(gen);
    }

    PhenotypeTransformer transformer;
    transformer.apply_iint(phenotype, covariates);

    auto mean = phenotype.mean();
    auto variance = (phenotype.array() - mean).square().mean();
    auto std_dev = std::sqrt(variance);

    REQUIRE_THAT(mean, Catch::Matchers::WithinAbs(0.0, 0.15));
    REQUIRE_THAT(std_dev, Catch::Matchers::WithinAbs(1.0, 0.15));
}

TEST_CASE(
    "PhenotypeTransformer - Custom offset parameter",
    "[phenotype_transformer]")
{
    Eigen::VectorXd phenotype(100);
    for (int i = 0; i < 100; ++i)
    {
        phenotype[i] = static_cast<double>(i);
    }

    auto phenotype_default = phenotype;
    auto phenotype_custom = phenotype;

    PhenotypeTransformer transformer_default(3.0 / 8.0);
    PhenotypeTransformer transformer_custom(0.5);

    transformer_default.apply_dint(phenotype_default);
    transformer_custom.apply_dint(phenotype_custom);

    REQUIRE((phenotype_default - phenotype_custom).norm() > 0.1);
}

TEST_CASE("PhenotypeTransformer - Small sample size", "[phenotype_transformer]")
{
    Eigen::VectorXd phenotype{{1.0, 2.0, 3.0, 4.0, 5.0}};

    PhenotypeTransformer transformer;

    REQUIRE_NOTHROW(transformer.apply_dint(phenotype));

    for (int i = 0; i < 4; ++i)
    {
        REQUIRE(phenotype[i] < phenotype[i + 1]);
    }
}

TEST_CASE("PhenotypeTransformer - Extreme values", "[phenotype_transformer]")
{
    Eigen::VectorXd phenotype(100);
    for (int i = 0; i < 100; ++i)
    {
        phenotype[i] = std::pow(10.0, static_cast<double>(i) / 10.0);
    }

    PhenotypeTransformer transformer;
    transformer.apply_dint(phenotype);

    REQUIRE(phenotype.allFinite());

    auto mean = phenotype.mean();
    auto std_dev = std::sqrt((phenotype.array() - mean).square().mean());

    REQUIRE_THAT(mean, Catch::Matchers::WithinAbs(0.0, 0.2));
    REQUIRE_THAT(std_dev, Catch::Matchers::WithinAbs(1.0, 0.2));
}
