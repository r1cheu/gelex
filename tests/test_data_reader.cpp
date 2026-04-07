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

#include <cstddef>

#include <catch2/catch_test_macros.hpp>

#include "file_fixture.h"
#include "gelex/data/dataframe/constants.h"
#include "gelex/data/reader.h"

using gelex::df::kSeparator;
using gelex::test::FileFixture;

TEST_CASE(
    "read_fam reads sample IDs as composite FID-IID index",
    "[data][reader]")
{
    FileFixture files;
    constexpr std::string_view kContent
        = "FAM1\tIID1\t0\t0\t1\t-9\n"
          "FAM1\tIID2\t0\t0\t2\t-9\n"
          "FAM2\tIID3\t0\t0\t1\t-9\n";

    auto path = files.create_text_file(kContent, ".fam");
    auto df = gelex::read_fam(path);

    REQUIRE(df.rows() == 3);
    REQUIRE(df.cols() == 0);

    auto key0 = std::string("FAM1") + kSeparator + "IID1";
    auto key1 = std::string("FAM1") + kSeparator + "IID2";
    auto key2 = std::string("FAM2") + kSeparator + "IID3";

    REQUIRE(df.index().at(key0) == 0);
    REQUIRE(df.index().at(key1) == 1);
    REQUIRE(df.index().at(key2) == 2);
}

TEST_CASE(
    "read_bim reads SNP ID as index with allele columns",
    "[data][reader]")
{
    FileFixture files;
    constexpr std::string_view kContent
        = "1\trs1\t0.0\t1000\tA\tG\n"
          "1\trs2\t0.0\t2000\tC\tT\n"
          "2\trs3\t0.0\t500\tA\tT\n";

    auto path = files.create_text_file(kContent, ".bim");
    auto df = gelex::read_bim(path);

    REQUIRE(df.rows() == 3);

    REQUIRE(df.index().at("rs1") == 0);
    REQUIRE(df.index().at("rs2") == 1);
    REQUIRE(df.index().at("rs3") == 2);

    auto a1 = df["A1"].as<std::string>();
    REQUIRE(a1[0] == "A");
    REQUIRE(a1[1] == "C");
    REQUIRE(a1[2] == "A");

    auto a2 = df["A2"].as<std::string>();
    REQUIRE(a2[0] == "G");
    REQUIRE(a2[1] == "T");
    REQUIRE(a2[2] == "T");
}

TEST_CASE(
    "read_pheno reads all phenotype columns with FID-IID index",
    "[data][reader]")
{
    FileFixture files;
    constexpr std::string_view kContent
        = "FID\tIID\tPheno1\tPheno2\n"
          "F1\tI1\t1.5\t2.5\n"
          "F1\tI2\t3.0\t4.0\n";

    auto path = files.create_text_file(kContent, ".tsv");
    auto df = gelex::read_pheno(path);

    REQUIRE(df.rows() == 2);
    REQUIRE(df.cols() == 2);
    REQUIRE(df["Pheno1"].as<double>()[0] == 1.5);
    REQUIRE(df["Pheno2"].as<double>()[1] == 4.0);
}

TEST_CASE(
    "read_pheno with pheno_col selects a single phenotype column",
    "[data][reader]")
{
    FileFixture files;
    constexpr std::string_view kContent
        = "FID\tIID\tPheno1\tPheno2\n"
          "F1\tI1\t1.5\t2.5\n"
          "F1\tI2\t3.0\t4.0\n";

    auto path = files.create_text_file(kContent, ".tsv");
    std::size_t col = 1;
    auto df = gelex::read_pheno(path, &col);

    REQUIRE(df.rows() == 2);
    REQUIRE(df.cols() == 1);
    REQUIRE(df["Pheno2"].as<double>()[0] == 2.5);
    REQUIRE(df["Pheno2"].as<double>()[1] == 4.0);
}

TEST_CASE(
    "read_qcovar reads quantitative covariates with FID-IID index",
    "[data][reader]")
{
    FileFixture files;
    constexpr std::string_view kContent
        = "FID\tIID\tAge\tHeight\n"
          "F1\tI1\t25.0\t170.5\n"
          "F1\tI2\t30.0\t165.0\n";

    auto path = files.create_text_file(kContent, ".tsv");
    auto df = gelex::read_qcovar(path);

    REQUIRE(df.rows() == 2);
    REQUIRE(df.cols() == 2);
    REQUIRE(df["Age"].as<double>()[0] == 25.0);
    REQUIRE(df["Height"].as<double>()[1] == 165.0);
}

TEST_CASE(
    "read_dcovar reads discrete covariates with FID-IID index",
    "[data][reader]")
{
    FileFixture files;
    constexpr std::string_view kContent
        = "FID\tIID\tSex\tRegion\n"
          "F1\tI1\tM\tNorth\n"
          "F1\tI2\tF\tSouth\n";

    auto path = files.create_text_file(kContent, ".tsv");
    auto df = gelex::read_dcovar(path);

    REQUIRE(df.rows() == 2);
    REQUIRE(df.cols() == 2);
    REQUIRE(df["Sex"].as<std::string>()[0] == "M");
    REQUIRE(df["Region"].as<std::string>()[1] == "South");
}
