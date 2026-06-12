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

#include <algorithm>
#include <limits>
#include <string>
#include <vector>

#include <Eigen/Core>
#include <catch2/catch_test_macros.hpp>

#include "bed_fixture.h"
#include "gelex/data/bed.h"
#include "gelex/data/dataframe/index.h"
#include "gelex/exception.h"

using gelex::GelexException;
using gelex::dataframe::Index;
using gelex::test::are_matrices_equal;
using gelex::test::BedFixture;

TEST_CASE("Bed reads BED SNP ranges", "[data][bed]")
{
    BedFixture fixture;
    const double nan = std::numeric_limits<double>::quiet_NaN();
    Eigen::MatrixXd genotypes{
        {0.0, 1.0, 2.0, nan},
        {2.0, 0.0, 1.0, 0.0},
        {1.0, 2.0, 0.0, 1.0},
        {0.0, 1.0, 1.0, 2.0},
    };

    const auto bed_prefix
        = fixture.create_deterministic_bed_files(genotypes).first;
    auto bed = gelex::open_bed(bed_prefix.string());

    REQUIRE(bed.num_samples() == genotypes.rows());
    REQUIRE(bed.num_snps() == genotypes.cols());
    REQUIRE(
        static_cast<Eigen::Index>(bed.sample_index().size())
        == genotypes.rows());
    REQUIRE(
        static_cast<Eigen::Index>(bed.snp_index().size()) == genotypes.cols());

    auto loaded = bed.read<double>();
    REQUIRE(are_matrices_equal(loaded, genotypes));

    Eigen::MatrixXd target_buf(genotypes.rows(), 2);
    bed.read_into<double>(target_buf, 1);

    Eigen::MatrixXd expected_chunk = genotypes.middleCols(1, 2);
    REQUIRE(are_matrices_equal(target_buf, expected_chunk));

    const auto float_chunk = bed.read<float>(0, 2);
    REQUIRE(float_chunk.cast<double>().isApprox(genotypes.leftCols(2)));
}

TEST_CASE("Bed reads into chunked output blocks", "[data][bed]")
{
    BedFixture fixture;
    Eigen::MatrixXd genotypes{
        {0.0, 1.0, 2.0, 0.0, 1.0},
        {2.0, 0.0, 1.0, 1.0, 2.0},
        {1.0, 2.0, 0.0, 2.0, 0.0},
        {0.0, 1.0, 1.0, 2.0, 1.0},
    };

    const auto bed_prefix
        = fixture.create_deterministic_bed_files(genotypes).first;
    auto bed = gelex::open_bed(bed_prefix.string());
    Eigen::MatrixXd loaded(genotypes.rows(), genotypes.cols());
    Eigen::MatrixXd buffer(genotypes.rows(), 2);

    for (Eigen::Index first = 0; first < bed.num_snps(); first += buffer.cols())
    {
        const auto width = std::min(buffer.cols(), bed.num_snps() - first);
        auto out = buffer.leftCols(width);

        bed.read_into<double>(out, first);
        loaded.middleCols(first, width) = out;
    }

    REQUIRE(loaded.isApprox(genotypes));
}

TEST_CASE("Bed reads selected SNP positions in target order", "[data][bed]")
{
    BedFixture fixture;
    Eigen::MatrixXd genotypes{
        {0.0, 1.0, 2.0, 0.0},
        {2.0, 0.0, 1.0, 1.0},
        {1.0, 2.0, 0.0, 2.0},
        {0.0, 1.0, 1.0, 2.0},
    };

    const auto bed_prefix
        = fixture.create_deterministic_bed_files(genotypes).first;
    auto bed = gelex::open_bed(bed_prefix.string());
    std::vector<Eigen::Index> selected{3, 1, 2};

    const auto loaded = bed.read_snps<double>(selected);
    Eigen::MatrixXd expected(genotypes.rows(), selected.size());
    expected.col(0) = genotypes.col(3);
    expected.col(1) = genotypes.col(1);
    expected.col(2) = genotypes.col(2);

    REQUIRE(loaded.isApprox(expected));

    Eigen::MatrixXd target_buf(genotypes.rows(), selected.size());
    bed.read_snps_into<double>(target_buf, selected);
    REQUIRE(target_buf.isApprox(expected));
}

TEST_CASE("Bed reads selected SNP IDs in target order", "[data][bed]")
{
    BedFixture fixture;
    Eigen::MatrixXd genotypes{
        {0.0, 1.0, 2.0, 0.0},
        {2.0, 0.0, 1.0, 1.0},
        {1.0, 2.0, 0.0, 2.0},
        {0.0, 1.0, 1.0, 2.0},
    };
    std::vector<std::string> snp_ids{"rs_a", "rs_b", "rs_c", "rs_d"};

    const auto bed_prefix
        = fixture.create_deterministic_bed_files(genotypes, {}, snp_ids).first;
    auto bed = gelex::open_bed(bed_prefix.string());
    Index<std::string> target_snps{
        std::vector<std::string>{"rs_d", "rs_b", "rs_c"}};

    const auto loaded = bed.read_snps<double>(target_snps);
    Eigen::MatrixXd expected(genotypes.rows(), 3);
    expected.col(0) = genotypes.col(3);
    expected.col(1) = genotypes.col(1);
    expected.col(2) = genotypes.col(2);

    REQUIRE(bed.snp_index().keys()[0] == snp_ids[0]);
    REQUIRE(bed.snp_index().keys()[1] == snp_ids[1]);
    REQUIRE(bed.snp_index().keys()[2] == snp_ids[2]);
    REQUIRE(bed.snp_index().keys()[3] == snp_ids[3]);
    REQUIRE(loaded.isApprox(expected));

    Eigen::MatrixXd target_buf(genotypes.rows(), target_snps.size());
    bed.read_snps_into<double>(target_buf, target_snps);
    REQUIRE(target_buf.isApprox(expected));
}

TEST_CASE("Bed projects samples to target index order", "[data][bed]")
{
    BedFixture fixture;
    Eigen::MatrixXd genotypes{
        {0.0, 1.0, 2.0, 0.0},
        {2.0, 0.0, 1.0, 1.0},
        {1.0, 2.0, 0.0, 2.0},
        {0.0, 1.0, 1.0, 2.0},
    };

    const auto bed_prefix
        = fixture.create_deterministic_bed_files(genotypes).first;
    auto source_bed = gelex::open_bed(bed_prefix.string());
    const auto source_keys = source_bed.sample_index().keys();
    std::vector<std::string> target_keys{
        source_keys[3],
        source_keys[1],
        source_keys[0],
    };
    Index<std::string> target_index{target_keys};

    auto projected_bed = gelex::open_bed(bed_prefix.string(), target_index);
    const auto loaded = projected_bed.read<double>(1, 4);
    const Eigen::MatrixXd expected{
        {genotypes(3, 1), genotypes(3, 2), genotypes(3, 3)},
        {genotypes(1, 1), genotypes(1, 2), genotypes(1, 3)},
        {genotypes(0, 1), genotypes(0, 2), genotypes(0, 3)},
    };

    REQUIRE(projected_bed.num_samples() == 3);
    REQUIRE(projected_bed.sample_index().keys()[0] == target_keys[0]);
    REQUIRE(projected_bed.sample_index().keys()[1] == target_keys[1]);
    REQUIRE(projected_bed.sample_index().keys()[2] == target_keys[2]);
    REQUIRE(loaded.isApprox(expected));
}

TEST_CASE("Bed rejects invalid ranges, SNPs, and target indices", "[data][bed]")
{
    BedFixture fixture;
    const auto bed_prefix = fixture.create_bed_files(4, 3, 0.0).first;
    auto bed = gelex::open_bed(bed_prefix.string());

    REQUIRE_THROWS_AS(bed.read<double>(-1, 1), GelexException);
    REQUIRE_THROWS_AS(bed.read<double>(2, 1), GelexException);
    REQUIRE_THROWS_AS(bed.read<double>(0, 4), GelexException);

    Eigen::MatrixXd wrong_rows(1, 1);
    REQUIRE_THROWS_AS(bed.read_into<double>(wrong_rows, 0), GelexException);

    Eigen::MatrixXd too_wide(bed.num_samples(), bed.num_snps() + 1);
    REQUIRE_THROWS_AS(bed.read_into<double>(too_wide, 0), GelexException);

    std::vector<Eigen::Index> bad_snps{0, bed.num_snps()};
    REQUIRE_THROWS_AS(bed.read_snps<double>(bad_snps), GelexException);

    Eigen::MatrixXd wrong_cols(bed.num_samples(), 1);
    std::vector<Eigen::Index> two_snps{0, 1};
    REQUIRE_THROWS_AS(
        bed.read_snps_into<double>(wrong_cols, two_snps), GelexException);

    Index<std::string> missing_snps{std::vector<std::string>{"missing_snp"}};
    REQUIRE_THROWS_AS(bed.read_snps<double>(missing_snps), GelexException);

    Index<std::string> missing_target{
        std::vector<std::string>{"missing_family_missing_sample"}};
    REQUIRE_THROWS_AS(
        gelex::open_bed(bed_prefix.string(), missing_target), GelexException);
}
