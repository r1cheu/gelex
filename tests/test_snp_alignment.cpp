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
#include <vector>

#include "gelex/data/bed.h"
#include "gelex/data/reader.h"
#include "gelex/data/snp_alignment.h"

#include "bed_fixture.h"
#include "file_fixture.h"

using gelex::test::BedFixture;
using gelex::test::FileFixture;

TEST_CASE("build_snp_alignment orients training SNPs onto bim", "[predict]")
{
    FileFixture files;

    SECTION("all SNPs match same orientation")
    {
        auto eff_path = files.create_text_file(
            "CHR\tSNP\tA1\tA2\tBETA_A\n"
            "1\trs1\tA\tG\t0.5\n"
            "1\trs2\tC\tT\t-0.2\n"
            "2\trs3\tA\tT\t0.8\n",
            ".snpeff");
        auto bim_path = files.create_text_file(
            "1\trs1\t0\t1000\tA\tG\n"
            "1\trs2\t0\t2000\tC\tT\n"
            "2\trs3\t0\t500\tA\tT\n",
            ".bim");

        auto plan = gelex::build_snp_alignment(
            gelex::read_snp_effects(eff_path), gelex::read_bim(bim_path));

        CHECK(plan.num_same == 3);
        CHECK(plan.num_flip == 0);
        CHECK(plan.missing_pos.empty());
        REQUIRE(plan.source_col.size() == 3);
        CHECK(plan.source_col == std::vector<Eigen::Index>{0, 1, 2});
        CHECK(plan.train_pos == std::vector<Eigen::Index>{0, 1, 2});
        CHECK(plan.flip == std::vector<char>{0, 0, 0});
    }

    SECTION("some SNPs absent from bim")
    {
        auto eff_path = files.create_text_file(
            "CHR\tSNP\tA1\tA2\tBETA_A\n"
            "1\trs1\tA\tG\t0.5\n"
            "1\trs2\tC\tT\t-0.2\n"
            "2\trs3\tA\tT\t0.8\n",
            ".snpeff");
        auto bim_path = files.create_text_file(
            "1\trs1\t0\t1000\tA\tG\n"
            "2\trs3\t0\t500\tA\tT\n",
            ".bim");

        auto plan = gelex::build_snp_alignment(
            gelex::read_snp_effects(eff_path), gelex::read_bim(bim_path));

        CHECK(plan.num_same == 2);
        CHECK(plan.num_absent == 1);
        CHECK(plan.source_col == std::vector<Eigen::Index>{0, 1});
        CHECK(plan.train_pos == std::vector<Eigen::Index>{0, 2});
        CHECK(plan.missing_pos == std::vector<Eigen::Index>{1});
    }

    SECTION(
        "swapped alleles orient as flip, unrelated alleles are incompatible")
    {
        auto eff_path = files.create_text_file(
            "CHR\tSNP\tA1\tA2\tBETA_A\n"
            "1\trs1\tA\tG\t0.5\n"
            "1\trs2\tC\tT\t-0.2\n"
            "2\trs3\tA\tT\t0.8\n",
            ".snpeff");
        auto bim_path = files.create_text_file(
            "1\trs1\t0\t1000\tA\tG\n"
            "1\trs2\t0\t2000\tT\tC\n"
            "2\trs3\t0\t500\tG\tT\n",
            ".bim");

        auto plan = gelex::build_snp_alignment(
            gelex::read_snp_effects(eff_path), gelex::read_bim(bim_path));

        CHECK(plan.num_same == 1);
        CHECK(plan.num_flip == 1);
        CHECK(plan.num_incompatible == 1);
        CHECK(plan.source_col == std::vector<Eigen::Index>{0, 1});
        CHECK(plan.train_pos == std::vector<Eigen::Index>{0, 1});
        CHECK(plan.flip == std::vector<char>{0, 1});
        CHECK(plan.missing_pos == std::vector<Eigen::Index>{2});
    }

    SECTION("no SNP matches")
    {
        auto eff_path = files.create_text_file(
            "CHR\tSNP\tA1\tA2\tBETA_A\n"
            "1\trs1\tA\tG\t0.5\n"
            "1\trs2\tC\tT\t-0.2\n",
            ".snpeff");
        auto bim_path = files.create_text_file(
            "1\trs4\t0\t1000\tA\tG\n"
            "1\trs5\t0\t2000\tC\tT\n",
            ".bim");

        auto plan = gelex::build_snp_alignment(
            gelex::read_snp_effects(eff_path), gelex::read_bim(bim_path));

        CHECK(plan.num_absent == 2);
        CHECK(plan.source_col.empty());
        CHECK(plan.train_pos.empty());
        CHECK(plan.missing_pos == std::vector<Eigen::Index>{0, 1});
    }
}

TEST_CASE(
    "load_aligned_genotypes orients and scatters onto the training axis",
    "[predict]")
{
    FileFixture files;
    BedFixture bed_files;

    // predict bfile: 2 samples x 3 SNPs, dosage counts bim A1
    Eigen::MatrixXd dosage{{2, 0, 1}, {1, 2, 0}};
    auto [prefix, _] = bed_files.create_deterministic_bed_files(
        dosage,
        {},
        {"rs1", "rs2", "rs3"},
        {},
        {{'A', 'G'}, {'C', 'T'}, {'A', 'T'}});

    // training axis: rs2 (same), rs4 (absent), rs1 (flipped G/A), rs3 (same)
    auto eff_path = files.create_text_file(
        "CHR\tSNP\tA1\tA2\tBETA_A\n"
        "1\trs2\tC\tT\t0.5\n"
        "1\trs4\tG\tA\t0.5\n"
        "1\trs1\tG\tA\t0.5\n"
        "2\trs3\tA\tT\t0.5\n",
        ".snpeff");

    auto plan = gelex::build_snp_alignment(
        gelex::read_snp_effects(eff_path),
        gelex::read_bim(prefix.string() + ".bim"));

    auto bed = gelex::open_bed(prefix.string());
    CHECK(plan.train_count == 4);
    auto genotype = gelex::load_aligned_genotypes(bed, plan);

    constexpr double NAN_VALUE = std::numeric_limits<double>::quiet_NaN();
    Eigen::MatrixXd expected{{0, NAN_VALUE, 0, 1}, {2, NAN_VALUE, 1, 0}};

    CHECK(gelex::test::are_matrices_equal(genotype, expected));
}
