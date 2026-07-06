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

#include <limits>
#include <optional>
#include <string>
#include <utility>
#include <vector>

#include <Eigen/Core>
#include <catch2/catch_test_macros.hpp>

#include "gelex/data/genotype_method.h"
#include "gelex/io/predict_reader.h"
#include "gelex/io/snpstats.h"
#include "gelex/predict/compute.h"
#include "gelex/predict/standardize.h"
#include "gelex/predict/types.h"

using gelex::Coefficients;
using gelex::GenotypeData;
using gelex::SnpEffects;
using gelex::SnpStats;
using gelex::SnpStatsData;
using gelex::standardize_genotypes;

namespace
{

auto make_snp_stats(
    gelex::GenotypeMethod method,
    Eigen::VectorXd mean,
    std::optional<Eigen::VectorXd> stddev = std::nullopt) -> SnpStats
{
    SnpStats stats;
    stats.method = method;
    stats.mean = std::move(mean);
    if (stddev)
    {
        stats.var = stddev->array().square();
    }
    return stats;
}

}  // namespace

// ---------------------------------------------------------------------------
// standardize_genotypes tests
// ---------------------------------------------------------------------------

TEST_CASE(
    "standardize_genotypes: additive StandardizeHWE",
    "[predict][standardize]")
{
    // 2x2 matrix: rows=samples, cols=SNPs
    // mean = [0.5, 1.5], stddev = [0.5, 0.5]
    // expected = (X - mean) / stddev
    GenotypeData geno;
    geno.add.resize(2, 2);
    geno.add << 0.0, 2.0, 1.0, 1.0;

    SnpStatsData snpstats;
    snpstats.add = make_snp_stats(
        gelex::GenotypeMethod::StandardizeHWE,
        Eigen::VectorXd{{0.5, 1.5}},
        Eigen::VectorXd{{0.5, 0.5}});
    snpstats.has_dom = false;

    standardize_genotypes(geno, snpstats);

    Eigen::MatrixXd expected(2, 2);
    expected << -1.0, 1.0, 1.0, -1.0;
    REQUIRE(geno.add.isApprox(expected));
}

TEST_CASE(
    "standardize_genotypes: additive CenterHWE (no stddev)",
    "[predict][standardize]")
{
    // mean-center only, no stddev division
    GenotypeData geno;
    geno.add.resize(2, 2);
    geno.add << 0.0, 2.0, 1.0, 1.0;

    SnpStatsData snpstats;
    snpstats.add = make_snp_stats(
        gelex::GenotypeMethod::CenterHWE, Eigen::VectorXd{{0.5, 1.5}});
    snpstats.has_dom = false;

    standardize_genotypes(geno, snpstats);

    Eigen::MatrixXd expected(2, 2);
    expected << -0.5, 0.5, 0.5, -0.5;
    REQUIRE(geno.add.isApprox(expected));
}

TEST_CASE(
    "standardize_genotypes: additive + dominance StandardizeHWE (RawPolicy)",
    "[predict][standardize]")
{
    // RawPolicy for Dom: 2 -> 0, others unchanged
    // add input: [[0, 2], [1, 1], [NaN, 0]]
    // add mean=[0.5, 1.0], stddev=[0.5, 1.0]
    // add expected: [[-1,1],[1,0],[0,-1]] (NaN -> 0 after center/scale)
    //
    // dom input same matrix; after RawPolicy: [[0, 0], [1, 1], [NaN, 0]]
    // dom mean=[0.5, 0.5], stddev=[0.5, 0.5]
    // dom expected: [[-1,-1],[1,1],[0,-1]]
    const double nan = std::numeric_limits<double>::quiet_NaN();

    GenotypeData geno;
    geno.add.resize(3, 2);
    geno.add << 0.0, 2.0, 1.0, 1.0, nan, 0.0;
    geno.dom = Eigen::MatrixXd(3, 2);
    (*geno.dom) << 0.0, 2.0, 1.0, 1.0, nan, 0.0;

    SnpStatsData snpstats;
    snpstats.add = make_snp_stats(
        gelex::GenotypeMethod::StandardizeHWE,
        Eigen::VectorXd{{0.5, 1.0}},
        Eigen::VectorXd{{0.5, 1.0}});
    snpstats.dom = make_snp_stats(
        gelex::GenotypeMethod::StandardizeHWE,
        Eigen::VectorXd{{0.5, 0.5}},
        Eigen::VectorXd{{0.5, 0.5}});
    snpstats.has_dom = true;

    standardize_genotypes(geno, snpstats);

    Eigen::MatrixXd expected_add(3, 2);
    expected_add << -1.0, 1.0, 1.0, 0.0, 0.0, -1.0;
    REQUIRE(geno.add.isApprox(expected_add));

    Eigen::MatrixXd expected_dom(3, 2);
    expected_dom << -1.0, -1.0, 1.0, 1.0, 0.0, -1.0;
    REQUIRE(geno.dom->isApprox(expected_dom));
}

TEST_CASE(
    "standardize_genotypes: additive + dominance OrthStandardizeHWE "
    "(OrthogonalPolicy)",
    "[predict][standardize]")
{
    // OrthogonalPolicy for Dom:
    // 2 -> 4*maf - 2
    // 1 -> 2*maf
    // 0 -> 0
    // maf = add_mean / 2
    //
    // add input: [[0, 2], [1, 1]]
    // add mean=[0.5, 1.5], stddev=[0.5, 0.5]
    // add expected: [[-1,1],[1,-1]]
    //
    // dom input: [[1, 0], [0, 1]]
    // maf0 = 0.5/2 = 0.25, maf1 = 1.5/2 = 0.75
    // OrthEncode col0: 1 -> 2*0.25=0.5, 0 -> 0.0 => [0.5, 0.0]
    // OrthEncode col1: 0 -> 0.0, 1 -> 2*0.75=1.5 => [0.0, 1.5]
    // dom mean=[0.25, 0.75], stddev=[0.25, 0.75]
    // dom expected col0: (0.5-0.25)/0.25=1.0, (0.0-0.25)/0.25=-1.0
    // dom expected col1: (0.0-0.75)/0.75=-1.0, (1.5-0.75)/0.75=1.0
    // dom expected: [[1,-1],[-1,1]]
    GenotypeData geno;
    geno.add.resize(2, 2);
    geno.add << 0.0, 2.0, 1.0, 1.0;
    geno.dom = Eigen::MatrixXd(2, 2);
    (*geno.dom) << 1.0, 0.0, 0.0, 1.0;

    SnpStatsData snpstats;
    snpstats.add = make_snp_stats(
        gelex::GenotypeMethod::OrthStandardizeHWE,
        Eigen::VectorXd{{0.5, 1.5}},
        Eigen::VectorXd{{0.5, 0.5}});
    snpstats.dom = make_snp_stats(
        gelex::GenotypeMethod::OrthStandardizeHWE,
        Eigen::VectorXd{{0.25, 0.75}},
        Eigen::VectorXd{{0.25, 0.75}});
    snpstats.has_dom = true;

    standardize_genotypes(geno, snpstats);

    Eigen::MatrixXd expected_add(2, 2);
    expected_add << -1.0, 1.0, 1.0, -1.0;
    REQUIRE(geno.add.isApprox(expected_add));

    Eigen::MatrixXd expected_dom(2, 2);
    expected_dom << 1.0, -1.0, -1.0, 1.0;
    REQUIRE(geno.dom->isApprox(expected_dom));
}

// ---------------------------------------------------------------------------
// compute_gebv tests
// ---------------------------------------------------------------------------

TEST_CASE("compute_gebv: additive only", "[predict][compute]")
{
    // add = [[1,2],[3,4]], effects.add = [0.5, -0.5]
    // add_predictions = [1*0.5 + 2*(-0.5), 3*0.5 + 4*(-0.5)] = [-0.5, -0.5]
    GenotypeData geno;
    geno.add.resize(2, 2);
    geno.add << 1.0, 2.0, 3.0, 4.0;

    SnpEffects effects;
    effects.add = Eigen::VectorXd{{0.5, -0.5}};

    const auto result = gelex::compute_gebv(geno, effects);

    Eigen::VectorXd expected_add{{-0.5, -0.5}};
    REQUIRE(result.add_predictions.isApprox(expected_add));
    REQUIRE(result.total.isApprox(expected_add));
    REQUIRE_FALSE(result.dom_predictions.has_value());
}

TEST_CASE("compute_gebv: additive + dominance", "[predict][compute]")
{
    // add = [[1,0],[0,1]], effects.add = [1.0, 2.0]
    // dom = [[0.5,0.5],[0.5,0.5]], effects.dom = [0.1, 0.2]
    // add_pred = [1*1.0+0*2.0, 0*1.0+1*2.0] = [1.0, 2.0]
    // dom_pred = [0.5*0.1+0.5*0.2, 0.5*0.1+0.5*0.2] = [0.15, 0.15]
    // total = [1.15, 2.15]
    GenotypeData geno;
    geno.add.resize(2, 2);
    geno.add << 1.0, 0.0, 0.0, 1.0;
    geno.dom = Eigen::MatrixXd(2, 2);
    (*geno.dom) << 0.5, 0.5, 0.5, 0.5;

    SnpEffects effects;
    effects.add = Eigen::VectorXd{{1.0, 2.0}};
    effects.dom = Eigen::VectorXd{{0.1, 0.2}};

    const auto result = gelex::compute_gebv(geno, effects);

    Eigen::VectorXd expected_add{{1.0, 2.0}};
    Eigen::VectorXd expected_dom{{0.15, 0.15}};
    Eigen::VectorXd expected_total{{1.15, 2.15}};
    REQUIRE(result.add_predictions.isApprox(expected_add));
    REQUIRE(result.dom_predictions->isApprox(expected_dom));
    REQUIRE(result.total.isApprox(expected_total));
}

// ---------------------------------------------------------------------------
// compute_covariate_effects tests
// ---------------------------------------------------------------------------

TEST_CASE("compute_covariate_effects: intercept only", "[predict][compute]")
{
    // covariates = [[1],[1],[1]], coefficients = {["Intercept"], [2.5]}
    // total = [2.5, 2.5, 2.5], per_covariate has 0 cols, covar_names empty
    Eigen::MatrixXd covariates(3, 1);
    covariates << 1.0, 1.0, 1.0;

    Coefficients coefficients;
    coefficients.names = {"Intercept"};
    coefficients.values = Eigen::VectorXd{{2.5}};

    const auto result
        = gelex::compute_covariate_effects(covariates, coefficients);

    Eigen::VectorXd expected_total{{2.5, 2.5, 2.5}};
    REQUIRE(result.total.isApprox(expected_total));
    REQUIRE(result.per_covariate.cols() == 0);
    REQUIRE(result.covar_names.empty());
}

TEST_CASE(
    "compute_covariate_effects: intercept + covariates",
    "[predict][compute]")
{
    // covariates = [[1, 25, 1], [1, 30, 0]]
    // coefficients = {["Intercept","Age","Sex_M"], [1.0, 0.2, -0.3]}
    // total = [1+5+(-0.3), 1+6+0] = [5.7, 7.0]
    // per_covariate: Age=[5.0, 6.0], Sex_M=[-0.3, 0.0]
    Eigen::MatrixXd covariates(2, 3);
    covariates << 1.0, 25.0, 1.0, 1.0, 30.0, 0.0;

    Coefficients coefficients;
    coefficients.names = {"Intercept", "Age", "Sex_M"};
    coefficients.values = Eigen::VectorXd{{1.0, 0.2, -0.3}};

    const auto result
        = gelex::compute_covariate_effects(covariates, coefficients);

    Eigen::VectorXd expected_total{{5.7, 7.0}};
    REQUIRE(result.total.isApprox(expected_total));

    REQUIRE(result.covar_names == std::vector<std::string>{"Age", "Sex_M"});
    REQUIRE(result.per_covariate.cols() == 2);

    Eigen::VectorXd expected_age{{5.0, 6.0}};
    Eigen::VectorXd expected_sex{{-0.3, 0.0}};
    REQUIRE(result.per_covariate.col(0).isApprox(expected_age));
    REQUIRE(result.per_covariate.col(1).isApprox(expected_sex));
}
