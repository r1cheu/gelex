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
        = "FAM1 IID1 0 0 1 -9\n"
          "FAM1 IID2 0 0 2 -9\n"
          "FAM2 IID3 0 0 1 -9\n";

    auto path = files.create_text_file(kContent, ".fam");
    auto df = gelex::read_fam(path);

    REQUIRE(df.rows() == 3);
    REQUIRE(df.cols() == 0);

    auto key0 = std::string("FAM1") + kSeparator + "IID1";
    auto key1 = std::string("FAM1") + kSeparator + "IID2";
    auto key2 = std::string("FAM2") + kSeparator + "IID3";

    REQUIRE(df.row_position(key0) == 0);
    REQUIRE(df.row_position(key1) == 1);
    REQUIRE(df.row_position(key2) == 2);
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

    REQUIRE(df.row_position("rs1") == 0);
    REQUIRE(df.row_position("rs2") == 1);
    REQUIRE(df.row_position("rs3") == 2);

    auto a1 = df.col("4").as<std::string>();
    REQUIRE(a1[0] == "A");
    REQUIRE(a1[1] == "C");
    REQUIRE(a1[2] == "A");

    auto a2 = df.col("5").as<std::string>();
    REQUIRE(a2[0] == "G");
    REQUIRE(a2[1] == "T");
    REQUIRE(a2[2] == "T");
}
