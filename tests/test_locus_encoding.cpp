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
#include <limits>
#include <vector>

#include <Eigen/Core>
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include "gelex/data/genotype_method.h"
#include "gelex/data/locus_encoding.h"
#include "gelex/data/locus_stats.h"
#include "gelex/types/genetic_mode.h"

using Catch::Matchers::WithinAbs;
using Catch::Matchers::WithinRel;
using namespace gelex;

namespace
{
constexpr double TOLERANCE = 1e-12;
}  // namespace

TEST_CASE("encoding_spec_from_method maps genotype methods", "[data][encoding]")
{
    SECTION("standardize uses empirical raw coding")
    {
        const gelex::EncodingSpec spec{gelex::encoding_spec_from_method(
            gelex::GeneticMode::D, gelex::GenotypeMethod::Standardize)};

        REQUIRE(spec.effect == gelex::GeneticMode::D);
        REQUIRE(spec.dominance_code == gelex::DominanceCode::Het);
        REQUIRE(spec.normalization == gelex::Normalization::CenterScale);
        REQUIRE(spec.moment_basis == gelex::MomentBasis::Empirical);
    }

    SECTION("center hwe uses theoretical raw coding")
    {
        const gelex::EncodingSpec spec{gelex::encoding_spec_from_method(
            gelex::GeneticMode::D, gelex::GenotypeMethod::CenterHWE)};

        REQUIRE(spec.effect == gelex::GeneticMode::D);
        REQUIRE(spec.dominance_code == gelex::DominanceCode::Het);
        REQUIRE(spec.normalization == gelex::Normalization::Center);
        REQUIRE(spec.moment_basis == gelex::MomentBasis::Theoretical);
    }

    SECTION("orthogonal standardize uses HWE dominance code")
    {
        const gelex::EncodingSpec spec{gelex::encoding_spec_from_method(
            gelex::GeneticMode::D, gelex::GenotypeMethod::OrthStandardize)};

        REQUIRE(spec.effect == gelex::GeneticMode::D);
        REQUIRE(spec.dominance_code == gelex::DominanceCode::HWE);
        REQUIRE(spec.normalization == gelex::Normalization::CenterScale);
        REQUIRE(spec.moment_basis == gelex::MomentBasis::Empirical);
    }

    SECTION("noia center uses NOIA dominance code")
    {
        const gelex::EncodingSpec spec{gelex::encoding_spec_from_method(
            gelex::GeneticMode::D, gelex::GenotypeMethod::NOIACenter)};

        REQUIRE(spec.effect == gelex::GeneticMode::D);
        REQUIRE(spec.dominance_code == gelex::DominanceCode::NOIA);
        REQUIRE(spec.normalization == gelex::Normalization::Center);
        REQUIRE(spec.moment_basis == gelex::MomentBasis::Empirical);
    }
}

TEST_CASE(
    "compute_locus_stats counts PLINK genotype classes",
    "[data][encoding]")
{
    const Eigen::VectorXd locus{
        {0.0, 1.0, 2.0, 1.0, std::numeric_limits<double>::quiet_NaN()}};

    const gelex::LocusStats stats{gelex::compute_locus_stats<double>(locus)};

    REQUIRE(stats.nA2A2 == 1);
    REQUIRE(stats.nA1A2 == 2);
    REQUIRE(stats.nA1A1 == 1);
    REQUIRE(stats.n_missing == 1);
    REQUIRE(stats.n_nonmissing() == 4);
    REQUIRE(stats.has_nonmissing());
    REQUIRE_THAT(stats.pA2A2(), WithinRel(0.25, TOLERANCE));
    REQUIRE_THAT(stats.pA1A2(), WithinRel(0.5, TOLERANCE));
    REQUIRE_THAT(stats.pA1A1(), WithinRel(0.25, TOLERANCE));
    REQUIRE_THAT(stats.A1freq(), WithinRel(0.5, TOLERANCE));
}

TEST_CASE(
    "make_locus_encoding fits additive center-scale code",
    "[data][encoding]")
{
    const gelex::LocusStats stats{
        .nA2A2 = 1, .nA1A2 = 2, .nA1A1 = 1, .n_missing = 1};
    const gelex::EncodingSpec spec{
        .effect = gelex::GeneticMode::A,
        .dominance_code = gelex::DominanceCode::NOIA,
        .normalization = gelex::Normalization::CenterScale,
        .moment_basis = gelex::MomentBasis::Empirical};

    const gelex::LocusEncoding encoding{
        gelex::detail::make_locus_encoding(7, stats, spec)};
    const double sd{std::sqrt(0.5)};

    REQUIRE(encoding.valid);
    REQUIRE(encoding.marker_index == 7);
    REQUIRE_THAT(encoding.mean, WithinRel(1.0, TOLERANCE));
    REQUIRE_THAT(encoding.var, WithinRel(0.5, TOLERANCE));
    REQUIRE_THAT(encoding.sd, WithinRel(sd, TOLERANCE));
    REQUIRE_THAT(encoding.code[0], WithinRel(-1.0 / sd, TOLERANCE));
    REQUIRE_THAT(encoding.code[1], WithinAbs(0.0, TOLERANCE));
    REQUIRE_THAT(encoding.code[2], WithinRel(1.0 / sd, TOLERANCE));
    REQUIRE_THAT(encoding.missing_encoded_value, WithinAbs(0.0, TOLERANCE));
}

TEST_CASE(
    "locus encoding additive modes match genotype processor values",
    "[data][encoding]")
{
    SECTION("empirical center-scale")
    {
        Eigen::MatrixXd genotypes{{0.0}, {1.0}, {2.0}, {1.0}, {0.0}};
        const EncodingSpec spec{encoding_spec_from_method(
            GeneticMode::A, GenotypeMethod::Standardize)};
        const LociEncoding encoding{
            detail::make_loci_encoding<double>(genotypes, spec)};
        const double sd{std::sqrt(0.56)};
        const Eigen::MatrixXd expected{
            {-0.8 / sd}, {0.2 / sd}, {1.2 / sd}, {0.2 / sd}, {-0.8 / sd}};

        detail::transform_inplace<double>(genotypes, encoding);

        REQUIRE(encoding.loci.size() == 1);
        REQUIRE(encoding.loci.front().valid);
        REQUIRE_THAT(encoding.loci.front().mean, WithinRel(0.8, TOLERANCE));
        REQUIRE_THAT(encoding.loci.front().sd, WithinRel(sd, TOLERANCE));
        REQUIRE(genotypes.isApprox(expected, TOLERANCE));
    }

    SECTION("empirical center")
    {
        Eigen::MatrixXd genotypes{{0.0}, {1.0}, {2.0}, {1.0}, {0.0}};
        const EncodingSpec spec{
            encoding_spec_from_method(GeneticMode::A, GenotypeMethod::Center)};
        const LociEncoding encoding{
            detail::make_loci_encoding<double>(genotypes, spec)};
        const Eigen::MatrixXd expected{{-0.8}, {0.2}, {1.2}, {0.2}, {-0.8}};

        detail::transform_inplace<double>(genotypes, encoding);

        REQUIRE(encoding.loci.front().valid);
        REQUIRE_THAT(encoding.loci.front().mean, WithinRel(0.8, TOLERANCE));
        REQUIRE(genotypes.isApprox(expected, TOLERANCE));
    }

    SECTION("theoretical center-scale")
    {
        Eigen::MatrixXd genotypes{{0.0}, {1.0}, {2.0}, {1.0}, {0.0}};
        const EncodingSpec spec{encoding_spec_from_method(
            GeneticMode::A, GenotypeMethod::StandardizeHWE)};
        const LociEncoding encoding{
            detail::make_loci_encoding<double>(genotypes, spec)};
        const double sd{std::sqrt(2.0 * 0.4 * 0.6)};
        const Eigen::MatrixXd expected{
            {-0.8 / sd}, {0.2 / sd}, {1.2 / sd}, {0.2 / sd}, {-0.8 / sd}};

        detail::transform_inplace<double>(genotypes, encoding);

        REQUIRE(encoding.loci.front().valid);
        REQUIRE_THAT(encoding.loci.front().mean, WithinRel(0.8, TOLERANCE));
        REQUIRE_THAT(encoding.loci.front().sd, WithinRel(sd, TOLERANCE));
        REQUIRE(genotypes.isApprox(expected, TOLERANCE));
    }
}

TEST_CASE(
    "locus encoding dominance modes match genotype processor values",
    "[data][encoding]")
{
    SECTION("heterozygote empirical center-scale")
    {
        Eigen::MatrixXd genotypes{{0.0}, {1.0}, {2.0}, {1.0}, {0.0}, {2.0}};
        const EncodingSpec spec{encoding_spec_from_method(
            GeneticMode::D, GenotypeMethod::Standardize)};
        const LociEncoding encoding{
            detail::make_loci_encoding<double>(genotypes, spec)};
        const double mean{1.0 / 3.0};
        const double sd{std::sqrt(2.0 / 9.0)};
        const Eigen::MatrixXd expected{
            {-mean / sd},
            {(1.0 - mean) / sd},
            {-mean / sd},
            {(1.0 - mean) / sd},
            {-mean / sd},
            {-mean / sd}};

        detail::transform_inplace<double>(genotypes, encoding);

        REQUIRE(encoding.loci.front().valid);
        REQUIRE_THAT(encoding.loci.front().mean, WithinRel(mean, TOLERANCE));
        REQUIRE_THAT(encoding.loci.front().sd, WithinRel(sd, TOLERANCE));
        REQUIRE(genotypes.isApprox(expected, TOLERANCE));
    }

    SECTION("heterozygote theoretical center-scale")
    {
        Eigen::MatrixXd genotypes{{0.0}, {1.0}, {2.0}, {1.0}, {0.0}};
        const EncodingSpec spec{encoding_spec_from_method(
            GeneticMode::D, GenotypeMethod::StandardizeHWE)};
        const LociEncoding encoding{
            detail::make_loci_encoding<double>(genotypes, spec)};
        const double mean{2.0 * 0.4 * 0.6};
        const double sd{
            std::sqrt(2.0 * 0.4 * 0.6 * ((0.4 * 0.4) + (0.6 * 0.6)))};
        const Eigen::MatrixXd expected{
            {-mean / sd},
            {(1.0 - mean) / sd},
            {-mean / sd},
            {(1.0 - mean) / sd},
            {-mean / sd}};

        detail::transform_inplace<double>(genotypes, encoding);

        REQUIRE(encoding.loci.front().valid);
        REQUIRE_THAT(encoding.loci.front().mean, WithinRel(mean, TOLERANCE));
        REQUIRE_THAT(encoding.loci.front().sd, WithinRel(sd, TOLERANCE));
        REQUIRE(genotypes.isApprox(expected, TOLERANCE));
    }

    SECTION("orthogonal empirical center")
    {
        Eigen::MatrixXd genotypes{{0.0}, {1.0}, {2.0}, {1.0}, {0.0}};
        const EncodingSpec spec{encoding_spec_from_method(
            GeneticMode::D, GenotypeMethod::OrthCenter)};
        const LociEncoding encoding{
            detail::make_loci_encoding<double>(genotypes, spec)};
        const Eigen::MatrixXd expected{
            {-0.24}, {0.56}, {-0.64}, {0.56}, {-0.24}};

        detail::transform_inplace<double>(genotypes, encoding);

        REQUIRE(encoding.loci.front().valid);
        REQUIRE_THAT(encoding.loci.front().mean, WithinRel(0.24, TOLERANCE));
        REQUIRE(genotypes.isApprox(expected, TOLERANCE));
    }

    SECTION("orthogonal empirical center-scale")
    {
        Eigen::MatrixXd genotypes{{0.0}, {1.0}, {2.0}, {1.0}, {0.0}};
        const EncodingSpec spec{encoding_spec_from_method(
            GeneticMode::D, GenotypeMethod::OrthStandardize)};
        const LociEncoding encoding{
            detail::make_loci_encoding<double>(genotypes, spec)};
        const double sd{std::sqrt(0.2304)};
        const Eigen::MatrixXd expected{
            {-0.24 / sd}, {0.56 / sd}, {-0.64 / sd}, {0.56 / sd}, {-0.24 / sd}};

        detail::transform_inplace<double>(genotypes, encoding);

        REQUIRE(encoding.loci.front().valid);
        REQUIRE_THAT(encoding.loci.front().mean, WithinRel(0.24, TOLERANCE));
        REQUIRE_THAT(encoding.loci.front().sd, WithinRel(sd, TOLERANCE));
        REQUIRE(genotypes.isApprox(expected, TOLERANCE));
    }

    SECTION("orthogonal theoretical center-scale")
    {
        Eigen::MatrixXd genotypes{{0.0}, {1.0}, {2.0}, {1.0}, {0.0}};
        const EncodingSpec spec{encoding_spec_from_method(
            GeneticMode::D, GenotypeMethod::OrthStandardizeHWE)};
        const LociEncoding encoding{
            detail::make_loci_encoding<double>(genotypes, spec)};
        const double sd{2.0 * 0.4 * 0.6};
        const Eigen::MatrixXd expected{
            {-0.32 / sd}, {0.48 / sd}, {-0.72 / sd}, {0.48 / sd}, {-0.32 / sd}};

        detail::transform_inplace<double>(genotypes, encoding);

        REQUIRE(encoding.loci.front().valid);
        REQUIRE_THAT(encoding.loci.front().mean, WithinRel(0.32, TOLERANCE));
        REQUIRE_THAT(encoding.loci.front().sd, WithinRel(sd, TOLERANCE));
        REQUIRE(genotypes.isApprox(expected, TOLERANCE));
    }
}

TEST_CASE(
    "locus encoding NOIA modes match genotype processor values",
    "[data][encoding]")
{
    SECTION("additive center-scale")
    {
        Eigen::MatrixXd genotypes{{0.0}, {1.0}, {2.0}, {1.0}, {0.0}};
        const EncodingSpec spec{encoding_spec_from_method(
            GeneticMode::A, GenotypeMethod::NOIAStandardize)};
        const LociEncoding encoding{
            detail::make_loci_encoding<double>(genotypes, spec)};
        const double sd{std::sqrt(0.56)};
        const Eigen::MatrixXd expected{
            {-0.8 / sd}, {0.2 / sd}, {1.2 / sd}, {0.2 / sd}, {-0.8 / sd}};

        detail::transform_inplace<double>(genotypes, encoding);

        REQUIRE(encoding.loci.front().valid);
        REQUIRE_THAT(encoding.loci.front().mean, WithinRel(0.8, TOLERANCE));
        REQUIRE_THAT(encoding.loci.front().sd, WithinRel(sd, TOLERANCE));
        REQUIRE(genotypes.isApprox(expected, TOLERANCE));
    }

    SECTION("dominance center")
    {
        Eigen::MatrixXd genotypes{{0.0}, {1.0}, {2.0}, {1.0}, {0.0}};
        const EncodingSpec spec{encoding_spec_from_method(
            GeneticMode::D, GenotypeMethod::NOIACenter)};
        const LociEncoding encoding{
            detail::make_loci_encoding<double>(genotypes, spec)};
        const double cA1A1{-0.32 / 0.56};
        const double cA1A2{0.32 / 0.56};
        const double cA2A2{-0.16 / 0.56};
        const Eigen::MatrixXd expected{
            {cA2A2}, {cA1A2}, {cA1A1}, {cA1A2}, {cA2A2}};

        detail::transform_inplace<double>(genotypes, encoding);

        REQUIRE(encoding.loci.front().valid);
        REQUIRE_THAT(encoding.loci.front().mean, WithinAbs(0.0, TOLERANCE));
        REQUIRE(genotypes.isApprox(expected, TOLERANCE));
    }

    SECTION("dominance center-scale")
    {
        Eigen::MatrixXd genotypes{{0.0}, {1.0}, {2.0}, {1.0}, {0.0}};
        const EncodingSpec spec{encoding_spec_from_method(
            GeneticMode::D, GenotypeMethod::NOIAStandardize)};
        const LociEncoding encoding{
            detail::make_loci_encoding<double>(genotypes, spec)};
        const double cA1A1{-0.32 / 0.56};
        const double cA1A2{0.32 / 0.56};
        const double cA2A2{-0.16 / 0.56};
        const double sd{std::sqrt(
            ((2.0 * cA2A2 * cA2A2) + (2.0 * cA1A2 * cA1A2) + (cA1A1 * cA1A1))
            / 5.0)};
        const Eigen::MatrixXd expected{
            {cA2A2 / sd},
            {cA1A2 / sd},
            {cA1A1 / sd},
            {cA1A2 / sd},
            {cA2A2 / sd}};

        detail::transform_inplace<double>(genotypes, encoding);

        REQUIRE(encoding.loci.front().valid);
        REQUIRE_THAT(encoding.loci.front().sd, WithinRel(sd, TOLERANCE));
        REQUIRE(genotypes.isApprox(expected, TOLERANCE));
    }

    SECTION("additive and dominance center codes are sample-orthogonal")
    {
        Eigen::MatrixXd additive{
            {2.0},
            {2.0},
            {2.0},
            {1.0},
            {0.0},
            {0.0},
            {0.0},
            {0.0},
            {0.0},
            {0.0}};
        Eigen::MatrixXd dominance{additive};
        const EncodingSpec additive_spec{encoding_spec_from_method(
            GeneticMode::A, GenotypeMethod::NOIACenter)};
        const EncodingSpec dominance_spec{encoding_spec_from_method(
            GeneticMode::D, GenotypeMethod::NOIACenter)};
        const LociEncoding additive_encoding{
            detail::make_loci_encoding<double>(additive, additive_spec)};
        const LociEncoding dominance_encoding{
            detail::make_loci_encoding<double>(dominance, dominance_spec)};

        detail::transform_inplace<double>(additive, additive_encoding);
        detail::transform_inplace<double>(dominance, dominance_encoding);

        REQUIRE_THAT(
            additive.col(0).dot(dominance.col(0)), WithinAbs(0.0, TOLERANCE));
    }
}

TEST_CASE(
    "locus encoding skips monomorphic and invalid loci",
    "[data][encoding]")
{
    SECTION("additive monomorphic")
    {
        Eigen::MatrixXd genotypes{{2.0}, {2.0}, {2.0}, {2.0}, {2.0}};
        const EncodingSpec spec{encoding_spec_from_method(
            GeneticMode::A, GenotypeMethod::Standardize)};
        const LociEncoding encoding{
            detail::make_loci_encoding<double>(genotypes, spec)};

        detail::transform_inplace<double>(genotypes, encoding);

        REQUIRE_FALSE(encoding.loci.front().valid);
        REQUIRE(genotypes.isZero(TOLERANCE));
    }

    SECTION("NOIA dominance without heterozygotes")
    {
        Eigen::MatrixXd genotypes{{0.0}, {2.0}, {0.0}, {2.0}};
        const EncodingSpec spec{encoding_spec_from_method(
            GeneticMode::D, GenotypeMethod::NOIAStandardize)};
        const LociEncoding encoding{
            detail::make_loci_encoding<double>(genotypes, spec)};

        detail::transform_inplace<double>(genotypes, encoding);

        REQUIRE_FALSE(encoding.loci.front().valid);
        REQUIRE(genotypes.isZero(TOLERANCE));
    }

    SECTION("missing genotype is encoded to analytic mean after centering")
    {
        const double nan_value{std::numeric_limits<double>::quiet_NaN()};
        Eigen::MatrixXd genotypes{{0.0}, {1.0}, {nan_value}, {2.0}, {1.0}};
        const EncodingSpec spec{encoding_spec_from_method(
            GeneticMode::A, GenotypeMethod::NOIAStandardize)};
        const LociEncoding encoding{
            detail::make_loci_encoding<double>(genotypes, spec)};

        detail::transform_inplace<double>(genotypes, encoding);

        REQUIRE(encoding.loci.front().valid);
        REQUIRE_THAT(genotypes(2, 0), WithinAbs(0.0, TOLERANCE));
    }
}

TEST_CASE(
    "make_loci_encoding keeps local marker indices without offset",
    "[data][encoding]")
{
    const Eigen::MatrixXd genotypes{
        {0.0, 2.0, 0.0},
        {1.0, 2.0, 0.0},
        {2.0, 2.0, 0.0},
        {1.0, 2.0, 0.0},
        {std::numeric_limits<double>::quiet_NaN(), 2.0, 0.0}};
    const gelex::EncodingSpec spec{
        .effect = gelex::GeneticMode::A,
        .dominance_code = gelex::DominanceCode::NOIA,
        .normalization = gelex::Normalization::CenterScale,
        .moment_basis = gelex::MomentBasis::Empirical};

    const gelex::LociEncoding encoding{
        gelex::detail::make_loci_encoding<double>(genotypes, spec)};

    REQUIRE(encoding.loci.size() == 3);
    REQUIRE(encoding.loci.front().column_index == 0);
    REQUIRE(encoding.loci.front().marker_index == 0);
    REQUIRE(encoding.loci.front().valid);
    REQUIRE(encoding.loci[1].column_index == 1);
    REQUIRE(encoding.loci[1].marker_index == 1);
    REQUIRE_FALSE(encoding.loci[1].valid);
    REQUIRE(encoding.loci[2].column_index == 2);
    REQUIRE(encoding.loci[2].marker_index == 2);
    REQUIRE_FALSE(encoding.loci[2].valid);
}

TEST_CASE(
    "make_loci_encoding separates chunk columns from marker indices",
    "[data][encoding]")
{
    Eigen::MatrixXd genotypes{
        {0.0, 2.0, 0.0},
        {1.0, 2.0, 0.0},
        {2.0, 2.0, 0.0},
        {1.0, 2.0, 0.0},
        {std::numeric_limits<double>::quiet_NaN(), 2.0, 0.0}};
    const gelex::EncodingSpec spec{
        .effect = gelex::GeneticMode::A,
        .dominance_code = gelex::DominanceCode::NOIA,
        .normalization = gelex::Normalization::CenterScale,
        .moment_basis = gelex::MomentBasis::Empirical};
    const gelex::LociEncoding encoding{
        gelex::detail::make_loci_encoding<double>(genotypes, spec, 1e-12, 10)};

    REQUIRE(encoding.loci.front().column_index == 0);
    REQUIRE(encoding.loci.front().marker_index == 10);
    REQUIRE(encoding.loci.front().valid);
    REQUIRE(encoding.loci[1].column_index == 1);
    REQUIRE(encoding.loci[1].marker_index == 11);
    REQUIRE_FALSE(encoding.loci[1].valid);
    REQUIRE(encoding.loci[2].column_index == 2);
    REQUIRE(encoding.loci[2].marker_index == 12);
    REQUIRE_FALSE(encoding.loci[2].valid);

    gelex::detail::transform_inplace<double>(genotypes, encoding);

    REQUIRE(genotypes.col(1).isZero(TOLERANCE));
    REQUIRE(genotypes.col(2).isZero(TOLERANCE));
}

TEST_CASE("transform_inplace zeroes skipped columns", "[data][encoding]")
{
    Eigen::MatrixXd genotypes{
        {0.0, 2.0, 0.0},
        {1.0, 2.0, 0.0},
        {2.0, 2.0, 0.0},
        {1.0, 2.0, 0.0},
        {std::numeric_limits<double>::quiet_NaN(), 2.0, 0.0}};
    const gelex::EncodingSpec spec{
        .effect = gelex::GeneticMode::A,
        .dominance_code = gelex::DominanceCode::NOIA,
        .normalization = gelex::Normalization::CenterScale,
        .moment_basis = gelex::MomentBasis::Empirical};
    const gelex::LociEncoding encoding{
        gelex::detail::make_loci_encoding<double>(genotypes, spec)};
    const double sd{std::sqrt(0.5)};
    const Eigen::MatrixXd expected{
        {-1.0 / sd, 0.0, 0.0},
        {0.0, 0.0, 0.0},
        {1.0 / sd, 0.0, 0.0},
        {0.0, 0.0, 0.0},
        {0.0, 0.0, 0.0}};

    gelex::detail::transform_inplace<double>(genotypes, encoding);

    REQUIRE(genotypes.isApprox(expected, TOLERANCE));
}

TEST_CASE("encode APIs fit and transform genotype matrices", "[data][encoding]")
{
    Eigen::MatrixXd inplace_genotypes{
        {0.0, 2.0},
        {1.0, 2.0},
        {2.0, 2.0},
        {1.0, 2.0},
        {std::numeric_limits<double>::quiet_NaN(), 2.0}};
    Eigen::MatrixXd raw{inplace_genotypes};
    Eigen::MatrixXd encoded{Eigen::MatrixXd::Zero(raw.rows(), raw.cols())};

    const gelex::LociEncoding inplace_encoding{gelex::encode_inplace<double>(
        inplace_genotypes,
        gelex::GeneticMode::A,
        gelex::GenotypeMethod::Standardize)};
    const gelex::LociEncoding into_encoding{gelex::encode_into<double, double>(
        raw,
        encoded,
        gelex::GeneticMode::A,
        gelex::GenotypeMethod::Standardize)};

    REQUIRE(inplace_encoding.loci.front().valid);
    REQUIRE_FALSE(inplace_encoding.loci[1].valid);
    REQUIRE(into_encoding.loci.front().valid);
    REQUIRE_FALSE(into_encoding.loci[1].valid);
    REQUIRE(inplace_genotypes.isApprox(encoded, TOLERANCE));
    REQUIRE(raw(0, 0) == 0.0);
    REQUIRE(raw(0, 1) == 2.0);
}
