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
#include <string>

#include <Eigen/Core>
#include <catch2/catch_test_macros.hpp>

#include "file_fixture.h"
#include "gelex/data/dataframe/constants.h"
#include "gelex/data/reader.h"
#include "gelex/io/predict/input_reader.h"
#include "gelex/predict/snp_alignment.h"

using gelex::dataframe::kSeparator;
using gelex::test::FileFixture;

TEST_CASE(
    "read_param reads term names and coefficients from TSV",
    "[predict][reader]")
{
    FileFixture files;
    constexpr std::string_view kContent
        = "term\tmean\n"
          "Intercept\t0.5\n"
          "Age\t-0.3\n"
          "Sex_M\t1.2\n";

    auto path = files.create_text_file(kContent, ".tsv");
    auto coefficients = gelex::predict::read_coefficients(path);

    REQUIRE(coefficients.names.size() == 3);
    REQUIRE(coefficients.names[0] == "Intercept");
    REQUIRE(coefficients.names[1] == "Age");
    REQUIRE(coefficients.names[2] == "Sex_M");

    REQUIRE(coefficients.values.size() == 3);
    REQUIRE(coefficients.values[0] == 0.5);
    REQUIRE(coefficients.values[1] == -0.3);
    REQUIRE(coefficients.values[2] == 1.2);
}

TEST_CASE(
    "read_covariates builds design matrix with intercept, qcovar, and dcovar",
    "[predict][reader]")
{
    FileFixture files;
    auto sep = std::string(1, kSeparator);

    // FAM file (space-delimited, no header): FID IID Father Mother Sex Pheno
    auto fam_path = files.create_text_file(
        "F1\tI1\t0\t0\t1\t-9\n"
        "F1\tI2\t0\t0\t2\t-9\n"
        "F2\tI1\t0\t0\t1\t-9\n",
        ".fam");

    // qcovar file (tab-delimited, header)
    auto qcovar_path = files.create_text_file(
        "FID\tIID\tAge\n"
        "F1\tI1\t25.0\n"
        "F1\tI2\t30.0\n"
        "F2\tI1\t35.0\n",
        ".tsv");

    // dcovar file (tab-delimited, header)
    auto dcovar_path = files.create_text_file(
        "FID\tIID\tSex\n"
        "F1\tI1\tM\n"
        "F1\tI2\tF\n"
        "F2\tI1\tM\n",
        ".tsv");

    auto sample_df = gelex::read_fam(fam_path);

    gelex::predict::Coefficients coefficients{
        .names = {"Intercept", "Age", "Sex" + sep + "M"},
        .values = Eigen::Vector3d{0.5, -0.3, 1.2},
    };

    auto covariates = gelex::predict::read_covariates(
        qcovar_path, dcovar_path, coefficients, sample_df);

    REQUIRE(covariates.rows() == 3);
    REQUIRE(covariates.cols() == 3);

    // col 0: Intercept (all 1s)
    REQUIRE(covariates(0, 0) == 1.0);
    REQUIRE(covariates(1, 0) == 1.0);
    REQUIRE(covariates(2, 0) == 1.0);

    // col 1: Age
    REQUIRE(covariates(0, 1) == 25.0);
    REQUIRE(covariates(1, 1) == 30.0);
    REQUIRE(covariates(2, 1) == 35.0);

    // col 2: Sex_M (M=1, F=0, M=1)
    REQUIRE(covariates(0, 2) == 1.0);
    REQUIRE(covariates(1, 2) == 0.0);
    REQUIRE(covariates(2, 2) == 1.0);
}

TEST_CASE(
    "read_snp_effects reads effect file with SNP ID as index",
    "[predict][reader]")
{
    FileFixture files;

    constexpr std::string_view kContent
        = "Chrom\tID\tPosition\tA1\tA2\tA1Freq\tAdd\tDom\n"
          "1\trs1\t1000\tA\tG\t0.3\t0.5\t0.1\n"
          "1\trs2\t2000\tC\tT\t0.4\t-0.2\t0.0\n"
          "2\trs3\t500\tA\tT\t0.1\t0.8\t-0.3\n";

    auto path = files.create_text_file(kContent, ".snp.eff");
    auto df = gelex::predict::read_snp_effects(path);

    REQUIRE(df.rows() == 3);

    REQUIRE(df.index().at("rs1") == 0);
    REQUIRE(df.index().at("rs2") == 1);
    REQUIRE(df.index().at("rs3") == 2);

    auto add = df["Add"].as<double>();
    REQUIRE(add[0] == 0.5);
    REQUIRE(add[1] == -0.2);
    REQUIRE(add[2] == 0.8);

    auto freq = df["A1Freq"].as<double>();
    REQUIRE(freq[0] == 0.3);
    REQUIRE(freq[1] == 0.4);
    REQUIRE(freq[2] == 0.1);

    REQUIRE(df.contains("Dom"));
    auto dom = df["Dom"].as<double>();
    REQUIRE(dom[0] == 0.1);
    REQUIRE(dom[1] == 0.0);
    REQUIRE(dom[2] == -0.3);
}

TEST_CASE(
    "build_snp_alignment maps effect SNPs to bim columns",
    "[predict][reader]")
{
    FileFixture files;

    SECTION("happy path: all SNPs match")
    {
        auto eff_path = files.create_text_file(
            "Chrom\tID\tA1\tA2\tAdd\n"
            "1\trs1\tA\tG\t0.5\n"
            "1\trs2\tC\tT\t-0.2\n"
            "2\trs3\tA\tT\t0.8\n",
            ".snp.eff");
        auto bim_path = files.create_text_file(
            "1\trs1\t0\t1000\tA\tG\n"
            "1\trs2\t0\t2000\tC\tT\n"
            "2\trs3\t0\t500\tA\tT\n",
            ".bim");

        auto snp_effects = gelex::predict::read_snp_effects(eff_path);
        auto bim_df = gelex::read_bim(bim_path);
        auto alignment
            = gelex::predict::build_snp_alignment(snp_effects, bim_df);

        REQUIRE(alignment.num_missing == 0);
        REQUIRE(alignment.num_mismatched == 0);
        REQUIRE(alignment.column_map.size() == 3);
        REQUIRE(alignment.column_map[0] == 0);
        REQUIRE(alignment.column_map[1] == 1);
        REQUIRE(alignment.column_map[2] == 2);
    }

    SECTION("partial match: some SNPs missing from bim")
    {
        auto eff_path = files.create_text_file(
            "Chrom\tID\tA1\tA2\tAdd\n"
            "1\trs1\tA\tG\t0.5\n"
            "1\trs2\tC\tT\t-0.2\n"
            "2\trs3\tA\tT\t0.8\n",
            ".snp.eff");
        auto bim_path = files.create_text_file(
            "1\trs1\t0\t1000\tA\tG\n"
            "2\trs3\t0\t500\tA\tT\n",
            ".bim");

        auto snp_effects = gelex::predict::read_snp_effects(eff_path);
        auto bim_df = gelex::read_bim(bim_path);
        auto alignment
            = gelex::predict::build_snp_alignment(snp_effects, bim_df);

        REQUIRE(alignment.num_missing == 1);
        REQUIRE(alignment.num_mismatched == 0);
        REQUIRE(alignment.column_map.size() == 3);
        REQUIRE(alignment.column_map[0] == 0);
        REQUIRE(alignment.column_map[1] == -1);
        REQUIRE(alignment.column_map[2] == 1);
    }

    SECTION("allele mismatch: ID matches but alleles differ")
    {
        auto eff_path = files.create_text_file(
            "Chrom\tID\tA1\tA2\tAdd\n"
            "1\trs1\tA\tG\t0.5\n"
            "1\trs2\tC\tT\t-0.2\n"
            "2\trs3\tA\tT\t0.8\n",
            ".snp.eff");
        auto bim_path = files.create_text_file(
            "1\trs1\t0\t1000\tA\tG\n"
            "1\trs2\t0\t2000\tT\tC\n"
            "2\trs3\t0\t500\tG\tT\n",
            ".bim");

        auto snp_effects = gelex::predict::read_snp_effects(eff_path);
        auto bim_df = gelex::read_bim(bim_path);
        auto alignment
            = gelex::predict::build_snp_alignment(snp_effects, bim_df);

        REQUIRE(alignment.num_missing == 0);
        REQUIRE(alignment.num_mismatched == 2);
        REQUIRE(alignment.column_map.size() == 3);
        REQUIRE(alignment.column_map[0] == 0);
        REQUIRE(alignment.column_map[1] == -1);
        REQUIRE(alignment.column_map[2] == -1);
    }

    SECTION("no match: all SNPs missing from bim")
    {
        auto eff_path = files.create_text_file(
            "Chrom\tID\tA1\tA2\tAdd\n"
            "1\trs1\tA\tG\t0.5\n"
            "1\trs2\tC\tT\t-0.2\n",
            ".snp.eff");
        auto bim_path = files.create_text_file(
            "1\trs4\t0\t1000\tA\tG\n"
            "1\trs5\t0\t2000\tC\tT\n",
            ".bim");

        auto snp_effects = gelex::predict::read_snp_effects(eff_path);
        auto bim_df = gelex::read_bim(bim_path);
        auto alignment
            = gelex::predict::build_snp_alignment(snp_effects, bim_df);

        REQUIRE(alignment.num_missing == 2);
        REQUIRE(alignment.num_mismatched == 0);
        REQUIRE(alignment.column_map.size() == 2);
        REQUIRE(alignment.column_map[0] == -1);
        REQUIRE(alignment.column_map[1] == -1);
    }
}
