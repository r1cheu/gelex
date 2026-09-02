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
#include <algorithm>
#include <catch2/catch_test_macros.hpp>
#include <string>
#include <vector>

#include "gelex/data/bed.h"
#include "gelex/data/dataframe/index.h"
#include "gelex/data/encode/encoder.h"
#include "gelex/data/encode/types.h"
#include "gelex/exception.h"
#include "gelex/genetic_mode.h"

#include "bed_fixture.h"

using gelex::DataFrameIndex;
using gelex::GelexException;
using gelex::test::BedFixture;

TEST_CASE("Bed exposes dataset metadata", "[data][bed]")
{
    BedFixture fixture;
    Eigen::MatrixXd genotypes{
        {0.0, 1.0, 2.0, 0.0},
        {2.0, 0.0, 1.0, 1.0},
        {1.0, 2.0, 0.0, 2.0},
        {0.0, 1.0, 1.0, 2.0},
    };
    std::vector<std::string> snp_ids{"rs_a", "rs_b", "rs_c", "rs_d"};

    const auto prefix
        = fixture.create_deterministic_bed_files(genotypes, {}, snp_ids).first;
    const auto bed = gelex::open_bed(prefix.string());

    REQUIRE(bed.num_samples() == genotypes.rows());
    REQUIRE(bed.num_snps() == genotypes.cols());
    REQUIRE(
        static_cast<Eigen::Index>(bed.sample_index().size())
        == genotypes.rows());
    REQUIRE(std::ranges::equal(bed.snp_index().keys(), snp_ids));
}

TEST_CASE("Bed gather projects samples into target order", "[data][bed]")
{
    BedFixture fixture;
    Eigen::MatrixXd genotypes{
        {0.0, 1.0, 2.0, 0.0},
        {2.0, 0.0, 1.0, 1.0},
        {1.0, 2.0, 0.0, 2.0},
        {0.0, 1.0, 1.0, 2.0},
    };

    const auto prefix = fixture.create_deterministic_bed_files(genotypes).first;
    auto bed = gelex::open_bed(prefix.string());

    const auto source_keys = bed.sample_index().keys();
    const std::vector<std::string> target_keys{
        source_keys[3], source_keys[1], source_keys[0]};
    bed.gather(DataFrameIndex<std::string>{target_keys});

    REQUIRE(bed.num_samples() == 3);
    REQUIRE(std::ranges::equal(bed.sample_index().keys(), target_keys));

    // The projection is exercised through its only consumer: expanding an
    // unnormalized additive encoding (value = dosage) straight from the packed
    // source must reproduce each SNP's genotypes in the projected target order.
    const gelex::LocusEncoder encoder{bed};
    const gelex::EncodingSpec spec{
        .effect = gelex::GeneticMode::A,
        .normalization = gelex::Normalization::None};

    for (Eigen::Index snp = 0; snp < bed.num_snps(); ++snp)
    {
        const auto stats = encoder.count(snp);
        const auto encoding = encoder.encoding(snp, stats, spec);

        Eigen::VectorXd got(bed.num_samples());
        encoder.expand(snp, encoding, got);

        const Eigen::VectorXd want{
            {genotypes(3, snp), genotypes(1, snp), genotypes(0, snp)}};
        INFO("snp=" << snp);
        CHECK(got.isApprox(want));
    }
}

TEST_CASE("Bed rejects an unknown gather target", "[data][bed]")
{
    BedFixture fixture;
    const auto prefix = fixture.create_bed_files(4, 3, 0.0).first;
    auto bed = gelex::open_bed(prefix.string());

    const DataFrameIndex<std::string> missing_target{
        std::vector<std::string>{"missing_family_missing_sample"}};
    REQUIRE_THROWS_AS(bed.gather(missing_target), GelexException);
}
