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
#include <string>
#include <string_view>

#include "gelex/data/reader.h"

#include "file_fixture.h"

using gelex::test::FileFixture;

TEST_CASE(
    "read_snp_effects reads effect file with SNP as index",
    "[data][reader]")
{
    FileFixture files;

    constexpr std::string_view CONTENT
        = "CHR\tSNP\tBP\tA1\tA2\tA1FREQ\tBETA_A\tBETA_D\n"
          "1\trs1\t1000\tA\tG\t0.3\t0.5\t0.1\n"
          "1\trs2\t2000\tC\tT\t0.4\t-0.2\t0.0\n"
          "2\trs3\t500\tA\tT\t0.1\t0.8\t-0.3\n";

    auto path = files.create_text_file(CONTENT, ".snpeff");
    auto df = gelex::read_snp_effects(path);

    REQUIRE(df.rows() == 3);

    REQUIRE(df.index().at("rs1") == 0);
    REQUIRE(df.index().at("rs2") == 1);
    REQUIRE(df.index().at("rs3") == 2);

    auto add = df["BETA_A"].as<double>();
    REQUIRE(add[0] == 0.5);
    REQUIRE(add[1] == -0.2);
    REQUIRE(add[2] == 0.8);

    auto freq = df["A1FREQ"].as<double>();
    REQUIRE(freq[0] == 0.3);
    REQUIRE(freq[1] == 0.4);
    REQUIRE(freq[2] == 0.1);

    REQUIRE(df.contains("BETA_D"));
    auto dom = df["BETA_D"].as<double>();
    REQUIRE(dom[0] == 0.1);
    REQUIRE(dom[1] == 0.0);
    REQUIRE(dom[2] == -0.3);
}
