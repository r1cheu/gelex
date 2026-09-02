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
#include <catch2/catch_test_macros.hpp>
#include <limits>

#include "gelex/data/bed.h"
#include "gelex/data/encode/encoder.h"
#include "gelex/data/encode/matrix.h"
#include "gelex/data/encode/spec.h"
#include "gelex/data/genotype_method.h"
#include "gelex/exception.h"
#include "gelex/genetic_mode.h"

#include "bed_fixture.h"

TEST_CASE(
    "Dense matrix and packed BED encoding remain equivalent",
    "[data][encode]")
{
    const double nan = std::numeric_limits<double>::quiet_NaN();
    const Eigen::MatrixXd raw{
        {0.0, 0.0, 2.0},
        {1.0, 0.0, 2.0},
        {2.0, 1.0, nan},
        {1.0, 2.0, 2.0},
        {nan, 2.0, 2.0}};
    gelex::test::BedFixture fixture;
    const auto files = fixture.create_deterministic_bed_files(raw);
    const gelex::Bed bed{gelex::open_bed(files.first.string())};
    const gelex::LocusEncoder encoder{bed};

    for (const gelex::GeneticMode mode :
         {gelex::GeneticMode::A, gelex::GeneticMode::D})
    {
        for (const auto [code, method] : gelex::GENOTYPE_METHOD_CODES)
        {
            Eigen::MatrixXd dense = raw;
            gelex::encode_inplace(dense, mode, method);

            Eigen::MatrixXd packed(raw.rows(), raw.cols());
            const gelex::EncodingSpec spec{
                gelex::encoding_spec_from_method(mode, method)};
            for (Eigen::Index marker = 0; marker < raw.cols(); ++marker)
            {
                const gelex::LocusStats stats{encoder.count(marker)};
                const gelex::LocusEncoding encoding{
                    encoder.encoding(marker, stats, spec)};
                encoder.expand(marker, encoding, packed.col(marker));
            }

            INFO("mode=" << static_cast<int>(mode) << " method=" << code);
            REQUIRE(dense.isApprox(packed));
        }
    }
}

TEST_CASE("Dense genotype matrices use Gelex locus encoding", "[data][encode]")
{
    using gelex::encode_inplace;
    using gelex::GeneticMode;
    using gelex::GenotypeMethod;

    SECTION("additive dosage is centered by empirical moments")
    {
        Eigen::MatrixXd genotypes{{0.0}, {1.0}, {2.0}};
        const Eigen::MatrixXd expected{{-1.0}, {0.0}, {1.0}};

        encode_inplace(genotypes, GeneticMode::A, GenotypeMethod::Center);

        REQUIRE(genotypes.isApprox(expected));
    }

    SECTION("dominance indicators are centered and missing values are zero")
    {
        const double nan = std::numeric_limits<double>::quiet_NaN();
        Eigen::MatrixXd genotypes{{0.0}, {1.0}, {2.0}, {1.0}, {nan}};
        const Eigen::MatrixXd expected{{-0.5}, {0.5}, {-0.5}, {0.5}, {0.0}};

        encode_inplace(genotypes, GeneticMode::D, GenotypeMethod::Center);

        REQUIRE(genotypes.isApprox(expected));
    }

    SECTION("invalid loci are zero")
    {
        Eigen::MatrixXd genotypes{{2.0, 0.0}, {2.0, 2.0}, {2.0, 0.0}};
        const Eigen::MatrixXd expected = Eigen::MatrixXd::Zero(3, 2);

        encode_inplace(
            genotypes, GeneticMode::D, GenotypeMethod::NOIAStandardize);

        REQUIRE(genotypes.isApprox(expected));
    }
}

TEST_CASE("Dense genotype matrices reject invalid dosage", "[data][encode]")
{
    Eigen::MatrixXd genotypes{{0.0}, {3.0}, {2.0}};

    REQUIRE_THROWS_AS(
        gelex::encode_inplace(
            genotypes,
            gelex::GeneticMode::A,
            gelex::GenotypeMethod::Standardize),
        gelex::GelexException);
}

TEST_CASE("Dense genotype encoding accepts empty matrices", "[data][encode]")
{
    Eigen::MatrixXd no_samples(0, 2);
    Eigen::MatrixXd no_markers(3, 0);

    gelex::encode_inplace(
        no_samples, gelex::GeneticMode::A, gelex::GenotypeMethod::Standardize);
    gelex::encode_inplace(
        no_markers, gelex::GeneticMode::A, gelex::GenotypeMethod::Standardize);

    REQUIRE(no_samples.rows() == 0);
    REQUIRE(no_samples.cols() == 2);
    REQUIRE(no_markers.rows() == 3);
    REQUIRE(no_markers.cols() == 0);
}
