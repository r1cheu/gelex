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

#include <vector>

#include <Eigen/Core>
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_exception.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>

#include "bed_fixture.h"
#include "file_fixture.h"

#include "gelex/data/reader/snp_effect_reader.h"
#include "gelex/exception.h"
#include "gelex/infra/logging/predict_event.h"
#include "gelex/predict/genotype_selector.h"

namespace fs = std::filesystem;

using namespace gelex;
using Catch::Matchers::ContainsSubstring;
using Catch::Matchers::WithinAbs;
using gelex::test::are_matrices_equal;
using gelex::test::FileFixture;

namespace
{

auto create_snp_effects(
    FileFixture& files,
    const std::string& header,
    const std::vector<std::string>& rows) -> SnpEffects
{
    std::string content = header + "\n";
    for (const auto& row : rows)
    {
        content += row + "\n";
    }
    auto file_path = files.create_text_file(content, ".snp.eff");
    gelex::detail::SnpEffectReader reader(file_path);
    return std::move(reader).take_effects();
}

}  // namespace

TEST_CASE("GenotypeSelector - perfect match", "[predict][genotype_selector]")
{
    FileFixture files;

    std::string bim_content
        = "1\trs001\t0\t1000\tA\tC\n"
          "1\trs002\t0\t2000\tT\tG\n"
          "1\trs003\t0\t3000\tC\tA\n";
    auto bim_path = files.create_text_file(bim_content, ".bim");
    auto bed_path = bim_path;
    bed_path.replace_extension(".bed");

    auto snp_effects = create_snp_effects(
        files,
        "Chrom\tPosition\tID\tA1\tA2\tA1Freq\tAdd\tDom",
        {"1\t1000\trs001\tA\tC\t0.25\t0.123\t0.045",
         "1\t2000\trs002\tT\tG\t0.75\t-0.456\t0.089",
         "1\t3000\trs003\tC\tA\t0.50\t0.789\t-0.012"});

    GenotypeSelector selector(bed_path, snp_effects);

    Eigen::MatrixXd genotypes(3, 3);
    genotypes << 0.0, 1.0, 2.0, 1.0, 2.0, 0.0, 2.0, 0.0, 1.0;

    Eigen::MatrixXd expected = genotypes;
    Eigen::MatrixXd result = selector.select(Eigen::MatrixXd(genotypes));

    REQUIRE(result.rows() == 3);
    REQUIRE(result.cols() == 3);
    REQUIRE(are_matrices_equal(result, expected, 1e-8));
}

TEST_CASE(
    "GenotypeSelector - column reordering",
    "[predict][genotype_selector]")
{
    FileFixture files;

    // BIM order: rs003, rs001, rs002
    std::string bim_content
        = "1\trs003\t0\t3000\tC\tA\n"
          "1\trs001\t0\t1000\tA\tC\n"
          "1\trs002\t0\t2000\tT\tG\n";
    auto bim_path = files.create_text_file(bim_content, ".bim");
    auto bed_path = bim_path;
    bed_path.replace_extension(".bed");

    // Effect order: rs001, rs002, rs003
    auto snp_effects = create_snp_effects(
        files,
        "Chrom\tPosition\tID\tA1\tA2\tA1Freq\tAdd\tDom",
        {"1\t1000\trs001\tA\tC\t0.25\t0.123\t0.045",
         "1\t2000\trs002\tT\tG\t0.75\t-0.456\t0.089",
         "1\t3000\trs003\tC\tA\t0.50\t0.789\t-0.012"});

    GenotypeSelector selector(bed_path, snp_effects);

    // Raw genotype columns match BIM order: rs003, rs001, rs002
    Eigen::MatrixXd genotypes(2, 3);
    genotypes << 0.5, 1.0, 1.5, 2.0, 0.0, 0.5;

    Eigen::MatrixXd result = selector.select(std::move(genotypes));

    // Result columns should match effect order: rs001, rs002, rs003
    // rs001 was BIM col 1, rs002 was BIM col 2, rs003 was BIM col 0
    REQUIRE(result.rows() == 2);
    REQUIRE(result.cols() == 3);
    REQUIRE(result(0, 0) == 1.0);  // rs001: BIM col 1
    REQUIRE(result(0, 1) == 1.5);  // rs002: BIM col 2
    REQUIRE(result(0, 2) == 0.5);  // rs003: BIM col 0
    REQUIRE(result(1, 0) == 0.0);
    REQUIRE(result(1, 1) == 0.5);
    REQUIRE(result(1, 2) == 2.0);
}

TEST_CASE(
    "GenotypeSelector - allele swap throws",
    "[predict][genotype_selector]")
{
    FileFixture files;

    // BIM has A1=C, A2=A (swapped relative to effect)
    std::string bim_content = "1\trs001\t0\t1000\tC\tA\n";
    auto bim_path = files.create_text_file(bim_content, ".bim");
    auto bed_path = bim_path;
    bed_path.replace_extension(".bed");

    auto snp_effects = create_snp_effects(
        files,
        "Chrom\tPosition\tID\tA1\tA2\tA1Freq\tAdd\tDom",
        {"1\t1000\trs001\tA\tC\t0.25\t0.123\t0.045"});

    REQUIRE_THROWS_AS(
        GenotypeSelector(bed_path, snp_effects), InvalidInputException);

    REQUIRE_THROWS_WITH(
        GenotypeSelector(bed_path, snp_effects),
        ContainsSubstring("Allele mismatch") && ContainsSubstring("rs001")
            && ContainsSubstring("plink2"));
}

TEST_CASE(
    "GenotypeSelector - incompatible alleles throw",
    "[predict][genotype_selector]")
{
    FileFixture files;

    // BIM has A1=A, A2=G (effect has A1=A, A2=C)
    std::string bim_content = "1\trs001\t0\t1000\tA\tG\n";
    auto bim_path = files.create_text_file(bim_content, ".bim");
    auto bed_path = bim_path;
    bed_path.replace_extension(".bed");

    auto snp_effects = create_snp_effects(
        files,
        "Chrom\tPosition\tID\tA1\tA2\tA1Freq\tAdd\tDom",
        {"1\t1000\trs001\tA\tC\t0.25\t0.123\t0.045"});

    REQUIRE_THROWS_AS(
        GenotypeSelector(bed_path, snp_effects), InvalidInputException);
}

TEST_CASE(
    "GenotypeSelector - few missing SNPs with mean imputation",
    "[predict][genotype_selector]")
{
    FileFixture files;

    // BIM has rs001 and rs003, missing rs002 (1/3 = 33%... need <=20%)
    // Use 10 SNPs with 2 missing = 20%
    std::string bim_content
        = "1\trs001\t0\t1000\tA\tC\n"
          "1\trs002\t0\t2000\tT\tG\n"
          "1\trs003\t0\t3000\tC\tA\n"
          "1\trs004\t0\t4000\tA\tG\n"
          "1\trs005\t0\t5000\tT\tC\n"
          "1\trs006\t0\t6000\tG\tA\n"
          "1\trs007\t0\t7000\tC\tT\n"
          "1\trs008\t0\t8000\tA\tT\n";
    auto bim_path = files.create_text_file(bim_content, ".bim");
    auto bed_path = bim_path;
    bed_path.replace_extension(".bed");

    // 10 effect SNPs, rs009 and rs010 are missing from BIM (2/10 = 20%)
    auto snp_effects = create_snp_effects(
        files,
        "Chrom\tPosition\tID\tA1\tA2\tA1Freq\tAdd\tDom",
        {"1\t1000\trs001\tA\tC\t0.25\t0.1\t0.01",
         "1\t2000\trs002\tT\tG\t0.75\t0.2\t0.02",
         "1\t3000\trs003\tC\tA\t0.50\t0.3\t0.03",
         "1\t4000\trs004\tA\tG\t0.40\t0.4\t0.04",
         "1\t5000\trs005\tT\tC\t0.60\t0.5\t0.05",
         "1\t6000\trs006\tG\tA\t0.30\t0.6\t0.06",
         "1\t7000\trs007\tC\tT\t0.20\t0.7\t0.07",
         "1\t8000\trs008\tA\tT\t0.10\t0.8\t0.08",
         "1\t9000\trs009\tG\tC\t0.35\t0.9\t0.09",
         "1\t10000\trs010\tT\tA\t0.45\t1.0\t0.10"});

    bool event_received = false;
    PredictObserver observer = [&](const PredictEvent& event)
    {
        std::visit(
            [&](const PredictSnpSelectionEvent& e)
            {
                event_received = true;
                REQUIRE(e.num_matched == 8);
                REQUIRE(e.num_missing == 2);
                REQUIRE(e.num_total == 10);
            },
            event);
    };

    GenotypeSelector selector(bed_path, snp_effects, observer);
    REQUIRE(event_received);

    // Raw genotype: 2 samples, 8 BIM columns
    Eigen::MatrixXd genotypes(2, 8);
    genotypes << 0.0, 1.0, 2.0, 1.0, 0.0, 2.0, 1.0, 0.0, 1.0, 0.0, 1.0, 2.0,
        1.0, 0.0, 2.0, 1.0;

    Eigen::MatrixXd result = selector.select(std::move(genotypes));

    REQUIRE(result.rows() == 2);
    REQUIRE(result.cols() == 10);

    // rs009 (index 8) missing, freq=0.35, mean = 2*0.35 = 0.70
    REQUIRE_THAT(result(0, 8), WithinAbs(0.70, 1e-8));
    REQUIRE_THAT(result(1, 8), WithinAbs(0.70, 1e-8));

    // rs010 (index 9) missing, freq=0.45, mean = 2*0.45 = 0.90
    REQUIRE_THAT(result(0, 9), WithinAbs(0.90, 1e-8));
    REQUIRE_THAT(result(1, 9), WithinAbs(0.90, 1e-8));

    // Matched SNPs should be copied correctly
    REQUIRE_THAT(result(0, 0), WithinAbs(0.0, 1e-8));  // rs001
    REQUIRE_THAT(result(0, 1), WithinAbs(1.0, 1e-8));  // rs002
}

TEST_CASE(
    "GenotypeSelector - too many missing SNPs throws",
    "[predict][genotype_selector]")
{
    FileFixture files;

    // BIM has only rs001
    std::string bim_content = "1\trs001\t0\t1000\tA\tC\n";
    auto bim_path = files.create_text_file(bim_content, ".bim");
    auto bed_path = bim_path;
    bed_path.replace_extension(".bed");

    // 4 effect SNPs, 3 missing = 75% > 20%
    auto snp_effects = create_snp_effects(
        files,
        "Chrom\tPosition\tID\tA1\tA2\tA1Freq\tAdd\tDom",
        {"1\t1000\trs001\tA\tC\t0.25\t0.123\t0.045",
         "1\t2000\trs002\tT\tG\t0.75\t-0.456\t0.089",
         "1\t3000\trs003\tC\tA\t0.50\t0.789\t-0.012",
         "1\t4000\trs004\tA\tG\t0.40\t0.111\t0.022"});

    REQUIRE_THROWS_AS(
        GenotypeSelector(bed_path, snp_effects), InvalidInputException);

    REQUIRE_THROWS_WITH(
        GenotypeSelector(bed_path, snp_effects),
        ContainsSubstring("Too many missing SNPs") && ContainsSubstring("3/4"));
}

TEST_CASE(
    "GenotypeSelector - case insensitive allele match",
    "[predict][genotype_selector]")
{
    FileFixture files;

    // BIM has lowercase alleles
    std::string bim_content
        = "1\trs001\t0\t1000\ta\tc\n"
          "1\trs002\t0\t2000\tt\tg\n";
    auto bim_path = files.create_text_file(bim_content, ".bim");
    auto bed_path = bim_path;
    bed_path.replace_extension(".bed");

    // Effect has uppercase alleles
    auto snp_effects = create_snp_effects(
        files,
        "Chrom\tPosition\tID\tA1\tA2\tA1Freq\tAdd\tDom",
        {"1\t1000\trs001\tA\tC\t0.25\t0.123\t0.045",
         "1\t2000\trs002\tT\tG\t0.75\t-0.456\t0.089"});

    REQUIRE_NOTHROW(GenotypeSelector(bed_path, snp_effects));

    GenotypeSelector selector(bed_path, snp_effects);

    Eigen::MatrixXd genotypes(2, 2);
    genotypes << 0.0, 1.0, 2.0, 0.0;

    Eigen::MatrixXd expected = genotypes;
    Eigen::MatrixXd result = selector.select(std::move(genotypes));

    REQUIRE(are_matrices_equal(result, expected, 1e-8));
}
