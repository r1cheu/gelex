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

#include <Eigen/Core>
#include <array>
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <cmath>
#include <cstddef>
#include <limits>

#include "gelex/data/encode/detail/encoding.h"
#include "gelex/data/encode/spec.h"
#include "gelex/data/encode/stats.h"
#include "gelex/data/encode/types.h"
#include "gelex/data/genotype_method.h"
#include "gelex/genetic_mode.h"

#include "locus_stats_oracle.h"

using Catch::Matchers::WithinAbs;
using Catch::Matchers::WithinRel;
using namespace gelex;

namespace
{
constexpr double tolerance = 1e-12;
constexpr std::array<Eigen::Index, 3> raw_code_by_dosage{3, 2, 0};

// Standardizes one genotype column through the fused per-locus path.
auto encode_column(const Eigen::VectorXd& column, const EncodingSpec& spec)
    -> LocusEncoding
{
    return detail::make_locus_encoding(
        0, gelex::test::compute_locus_stats<double>(column), spec);
}
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

    const gelex::LocusStats stats{
        gelex::test::compute_locus_stats<double>(locus)};

    REQUIRE(stats.nA2A2 == 1);
    REQUIRE(stats.nA1A2 == 2);
    REQUIRE(stats.nA1A1 == 1);
    REQUIRE(stats.n_missing == 1);
    REQUIRE(stats.n_nonmissing() == 4);
    REQUIRE(stats.has_nonmissing());
    REQUIRE_THAT(stats.pA2A2(), WithinRel(0.25, tolerance));
    REQUIRE_THAT(stats.pA1A2(), WithinRel(0.5, tolerance));
    REQUIRE_THAT(stats.pA1A1(), WithinRel(0.25, tolerance));
    REQUIRE_THAT(stats.A1freq(), WithinRel(0.5, tolerance));
}

TEST_CASE(
    "make_locus_encoding fits additive center-scale LUT",
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
    REQUIRE_THAT(encoding.mean, WithinRel(1.0, tolerance));
    REQUIRE_THAT(encoding.var, WithinRel(0.5, tolerance));
    REQUIRE_THAT(encoding.sd, WithinRel(sd, tolerance));
    REQUIRE_THAT(encoding.lut[3], WithinRel(-1.0 / sd, tolerance));
    REQUIRE_THAT(encoding.lut[2], WithinAbs(0.0, tolerance));
    REQUIRE_THAT(encoding.lut[0], WithinRel(1.0 / sd, tolerance));
    REQUIRE_THAT(encoding.lut[1], WithinAbs(0.0, tolerance));
}

TEST_CASE(
    "additive locus LUTs match genotype processor values",
    "[data][encoding]")
{
    const Eigen::VectorXd genotypes{{0.0, 1.0, 2.0, 1.0, 0.0}};

    SECTION("empirical center-scale")
    {
        const EncodingSpec spec{encoding_spec_from_method(
            GeneticMode::A, GenotypeMethod::Standardize)};
        const LocusEncoding enc{encode_column(genotypes, spec)};
        const double sd{std::sqrt(0.56)};

        REQUIRE(enc.valid);
        REQUIRE_THAT(enc.mean, WithinRel(0.8, tolerance));
        REQUIRE_THAT(enc.sd, WithinRel(sd, tolerance));
        REQUIRE_THAT(enc.lut[3], WithinRel(-0.8 / sd, tolerance));
        REQUIRE_THAT(enc.lut[2], WithinRel(0.2 / sd, tolerance));
        REQUIRE_THAT(enc.lut[0], WithinRel(1.2 / sd, tolerance));
    }

    SECTION("empirical center")
    {
        const EncodingSpec spec{
            encoding_spec_from_method(GeneticMode::A, GenotypeMethod::Center)};
        const LocusEncoding enc{encode_column(genotypes, spec)};

        REQUIRE(enc.valid);
        REQUIRE_THAT(enc.mean, WithinRel(0.8, tolerance));
        REQUIRE_THAT(enc.lut[3], WithinRel(-0.8, tolerance));
        REQUIRE_THAT(enc.lut[2], WithinRel(0.2, tolerance));
        REQUIRE_THAT(enc.lut[0], WithinRel(1.2, tolerance));
    }

    SECTION("theoretical center-scale")
    {
        const EncodingSpec spec{encoding_spec_from_method(
            GeneticMode::A, GenotypeMethod::StandardizeHWE)};
        const LocusEncoding enc{encode_column(genotypes, spec)};
        const double sd{std::sqrt(2.0 * 0.4 * 0.6)};

        REQUIRE(enc.valid);
        REQUIRE_THAT(enc.mean, WithinRel(0.8, tolerance));
        REQUIRE_THAT(enc.sd, WithinRel(sd, tolerance));
        REQUIRE_THAT(enc.lut[3], WithinRel(-0.8 / sd, tolerance));
        REQUIRE_THAT(enc.lut[2], WithinRel(0.2 / sd, tolerance));
        REQUIRE_THAT(enc.lut[0], WithinRel(1.2 / sd, tolerance));
    }
}

TEST_CASE(
    "dominance locus LUTs match genotype processor values",
    "[data][encoding]")
{
    SECTION("heterozygote empirical center-scale")
    {
        const Eigen::VectorXd genotypes{{0.0, 1.0, 2.0, 1.0, 0.0, 2.0}};
        const EncodingSpec spec{encoding_spec_from_method(
            GeneticMode::D, GenotypeMethod::Standardize)};
        const LocusEncoding enc{encode_column(genotypes, spec)};
        const double mean{1.0 / 3.0};
        const double sd{std::sqrt(2.0 / 9.0)};

        REQUIRE(enc.valid);
        REQUIRE_THAT(enc.mean, WithinRel(mean, tolerance));
        REQUIRE_THAT(enc.sd, WithinRel(sd, tolerance));
        REQUIRE_THAT(enc.lut[3], WithinRel(-mean / sd, tolerance));
        REQUIRE_THAT(enc.lut[2], WithinRel((1.0 - mean) / sd, tolerance));
        REQUIRE_THAT(enc.lut[0], WithinRel(-mean / sd, tolerance));
    }

    SECTION("heterozygote theoretical center-scale")
    {
        const Eigen::VectorXd genotypes{{0.0, 1.0, 2.0, 1.0, 0.0}};
        const EncodingSpec spec{encoding_spec_from_method(
            GeneticMode::D, GenotypeMethod::StandardizeHWE)};
        const LocusEncoding enc{encode_column(genotypes, spec)};
        const double mean{2.0 * 0.4 * 0.6};
        const double sd{
            std::sqrt(2.0 * 0.4 * 0.6 * ((0.4 * 0.4) + (0.6 * 0.6)))};

        REQUIRE(enc.valid);
        REQUIRE_THAT(enc.mean, WithinRel(mean, tolerance));
        REQUIRE_THAT(enc.sd, WithinRel(sd, tolerance));
        REQUIRE_THAT(enc.lut[3], WithinRel(-mean / sd, tolerance));
        REQUIRE_THAT(enc.lut[2], WithinRel((1.0 - mean) / sd, tolerance));
        REQUIRE_THAT(enc.lut[0], WithinRel(-mean / sd, tolerance));
    }

    SECTION("orthogonal empirical center")
    {
        const Eigen::VectorXd genotypes{{0.0, 1.0, 2.0, 1.0, 0.0}};
        const EncodingSpec spec{encoding_spec_from_method(
            GeneticMode::D, GenotypeMethod::OrthCenter)};
        const LocusEncoding enc{encode_column(genotypes, spec)};

        REQUIRE(enc.valid);
        REQUIRE_THAT(enc.mean, WithinRel(0.24, tolerance));
        REQUIRE_THAT(enc.lut[3], WithinRel(-0.24, tolerance));
        REQUIRE_THAT(enc.lut[2], WithinRel(0.56, tolerance));
        REQUIRE_THAT(enc.lut[0], WithinRel(-0.64, tolerance));
    }

    SECTION("orthogonal empirical center-scale")
    {
        const Eigen::VectorXd genotypes{{0.0, 1.0, 2.0, 1.0, 0.0}};
        const EncodingSpec spec{encoding_spec_from_method(
            GeneticMode::D, GenotypeMethod::OrthStandardize)};
        const LocusEncoding enc{encode_column(genotypes, spec)};
        const double sd{std::sqrt(0.2304)};

        REQUIRE(enc.valid);
        REQUIRE_THAT(enc.mean, WithinRel(0.24, tolerance));
        REQUIRE_THAT(enc.sd, WithinRel(sd, tolerance));
        REQUIRE_THAT(enc.lut[3], WithinRel(-0.24 / sd, tolerance));
        REQUIRE_THAT(enc.lut[2], WithinRel(0.56 / sd, tolerance));
        REQUIRE_THAT(enc.lut[0], WithinRel(-0.64 / sd, tolerance));
    }

    SECTION("orthogonal theoretical center-scale")
    {
        const Eigen::VectorXd genotypes{{0.0, 1.0, 2.0, 1.0, 0.0}};
        const EncodingSpec spec{encoding_spec_from_method(
            GeneticMode::D, GenotypeMethod::OrthStandardizeHWE)};
        const LocusEncoding enc{encode_column(genotypes, spec)};
        const double sd{2.0 * 0.4 * 0.6};

        REQUIRE(enc.valid);
        REQUIRE_THAT(enc.mean, WithinRel(0.32, tolerance));
        REQUIRE_THAT(enc.sd, WithinRel(sd, tolerance));
        REQUIRE_THAT(enc.lut[3], WithinRel(-0.32 / sd, tolerance));
        REQUIRE_THAT(enc.lut[2], WithinRel(0.48 / sd, tolerance));
        REQUIRE_THAT(enc.lut[0], WithinRel(-0.72 / sd, tolerance));
    }
}

TEST_CASE("NOIA locus LUTs match genotype processor values", "[data][encoding]")
{
    const double cA1A1{-0.32 / 0.56};
    const double cA1A2{0.32 / 0.56};
    const double cA2A2{-0.16 / 0.56};

    SECTION("additive center-scale")
    {
        const Eigen::VectorXd genotypes{{0.0, 1.0, 2.0, 1.0, 0.0}};
        const EncodingSpec spec{encoding_spec_from_method(
            GeneticMode::A, GenotypeMethod::NOIAStandardize)};
        const LocusEncoding enc{encode_column(genotypes, spec)};
        const double sd{std::sqrt(0.56)};

        REQUIRE(enc.valid);
        REQUIRE_THAT(enc.mean, WithinRel(0.8, tolerance));
        REQUIRE_THAT(enc.sd, WithinRel(sd, tolerance));
        REQUIRE_THAT(enc.lut[3], WithinRel(-0.8 / sd, tolerance));
        REQUIRE_THAT(enc.lut[2], WithinRel(0.2 / sd, tolerance));
        REQUIRE_THAT(enc.lut[0], WithinRel(1.2 / sd, tolerance));
    }

    SECTION("dominance center")
    {
        const Eigen::VectorXd genotypes{{0.0, 1.0, 2.0, 1.0, 0.0}};
        const EncodingSpec spec{encoding_spec_from_method(
            GeneticMode::D, GenotypeMethod::NOIACenter)};
        const LocusEncoding enc{encode_column(genotypes, spec)};

        REQUIRE(enc.valid);
        REQUIRE_THAT(enc.mean, WithinAbs(0.0, tolerance));
        REQUIRE_THAT(enc.lut[3], WithinRel(cA2A2, tolerance));
        REQUIRE_THAT(enc.lut[2], WithinRel(cA1A2, tolerance));
        REQUIRE_THAT(enc.lut[0], WithinRel(cA1A1, tolerance));
    }

    SECTION("dominance center-scale")
    {
        const Eigen::VectorXd genotypes{{0.0, 1.0, 2.0, 1.0, 0.0}};
        const EncodingSpec spec{encoding_spec_from_method(
            GeneticMode::D, GenotypeMethod::NOIAStandardize)};
        const LocusEncoding enc{encode_column(genotypes, spec)};
        const double sd{std::sqrt(
            ((2.0 * cA2A2 * cA2A2) + (2.0 * cA1A2 * cA1A2) + (cA1A1 * cA1A1))
            / 5.0)};

        REQUIRE(enc.valid);
        REQUIRE_THAT(enc.sd, WithinRel(sd, tolerance));
        REQUIRE_THAT(enc.lut[3], WithinRel(cA2A2 / sd, tolerance));
        REQUIRE_THAT(enc.lut[2], WithinRel(cA1A2 / sd, tolerance));
        REQUIRE_THAT(enc.lut[0], WithinRel(cA1A1 / sd, tolerance));
    }

    SECTION("additive and dominance center LUTs are sample-orthogonal")
    {
        const Eigen::VectorXd genotypes{
            {2.0, 2.0, 2.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}};
        const EncodingSpec additive_spec{encoding_spec_from_method(
            GeneticMode::A, GenotypeMethod::NOIACenter)};
        const EncodingSpec dominance_spec{encoding_spec_from_method(
            GeneticMode::D, GenotypeMethod::NOIACenter)};
        const LocusEncoding additive{encode_column(genotypes, additive_spec)};
        const LocusEncoding dominance{encode_column(genotypes, dominance_spec)};

        const auto apply = [&genotypes](const LocusEncoding& enc)
        {
            Eigen::VectorXd out(genotypes.size());
            for (Eigen::Index i = 0; i < genotypes.size(); ++i)
            {
                out[i] = enc.lut[raw_code_by_dosage[static_cast<std::size_t>(
                    genotypes[i])]];
            }
            return out;
        };

        REQUIRE_THAT(
            apply(additive).dot(apply(dominance)), WithinAbs(0.0, tolerance));
    }
}

TEST_CASE(
    "locus encoding flags monomorphic and invalid loci",
    "[data][encoding]")
{
    SECTION("additive monomorphic")
    {
        const Eigen::VectorXd genotypes{{2.0, 2.0, 2.0, 2.0, 2.0}};
        const EncodingSpec spec{encoding_spec_from_method(
            GeneticMode::A, GenotypeMethod::Standardize)};

        const LocusEncoding enc{encode_column(genotypes, spec)};
        REQUIRE_FALSE(enc.valid);
        REQUIRE(enc.lut.isZero());
    }

    SECTION("NOIA dominance without heterozygotes")
    {
        const Eigen::VectorXd genotypes{{0.0, 2.0, 0.0, 2.0}};
        const EncodingSpec spec{encoding_spec_from_method(
            GeneticMode::D, GenotypeMethod::NOIAStandardize)};

        const LocusEncoding enc{encode_column(genotypes, spec)};
        REQUIRE_FALSE(enc.valid);
        REQUIRE(enc.lut.isZero());
    }

    SECTION("centered encoding maps missing genotype to zero")
    {
        const double nan_value{std::numeric_limits<double>::quiet_NaN()};
        const Eigen::VectorXd genotypes{{0.0, 1.0, nan_value, 2.0, 1.0}};
        const EncodingSpec spec{encoding_spec_from_method(
            GeneticMode::A, GenotypeMethod::NOIAStandardize)};
        const LocusEncoding enc{encode_column(genotypes, spec)};

        REQUIRE(enc.valid);
        REQUIRE_THAT(enc.lut[1], WithinAbs(0.0, tolerance));
    }

    SECTION("unnormalized encoding maps missing genotype to the analytic mean")
    {
        const Eigen::VectorXd genotypes{
            {0.0, 1.0, std::numeric_limits<double>::quiet_NaN(), 2.0, 1.0}};
        const EncodingSpec spec{
            .effect = GeneticMode::A,
            .normalization = Normalization::None,
            .moment_basis = MomentBasis::Empirical};
        const LocusEncoding enc{encode_column(genotypes, spec)};

        REQUIRE(enc.valid);
        REQUIRE_THAT(enc.mean, WithinAbs(1.0, tolerance));
        REQUIRE_THAT(enc.lut[1], WithinAbs(enc.mean, tolerance));
    }
}
