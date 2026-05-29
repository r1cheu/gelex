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

#include <cstdint>
#include <string>

#include <catch2/catch_test_macros.hpp>

#include "file_fixture.h"
#include "gelex/data/reader.h"
#include "gelex/exception.h"

using gelex::test::FileFixture;

TEST_CASE("read_bim - Valid file parsing", "[reader][bim]")
{
    FileFixture files;

    SECTION("Happy path - tab-delimited file")
    {
        auto file_path = files.create_text_file(
            "1\trs12345\t0\t1000\tA\tG\n"
            "1\trs67890\t0.001\t2000\tC\tT\n"
            "2\trs24680\t0.002\t3000\tG\tA\n"
            "X\trs13579\t0.003\t4000\tT\tC\n"
            "1\trs11223\t0.004\t5000\tA\tT",
            ".bim");

        auto df = gelex::read_bim(file_path);

        REQUIRE(df.rows() == 5);

        auto keys = df.index().keys();
        REQUIRE(keys[0] == "rs12345");
        REQUIRE(keys[1] == "rs67890");
        REQUIRE(keys[2] == "rs24680");
        REQUIRE(keys[3] == "rs13579");
        REQUIRE(keys[4] == "rs11223");

        auto chrom = df["chrom"].as<std::string>();
        auto pos = df["pos"].as<std::int32_t>();
        auto a1 = df["A1"].as<std::string>();
        auto a2 = df["A2"].as<std::string>();

        REQUIRE(chrom[0] == "1");
        REQUIRE(pos[0] == 1000);
        REQUIRE(a1[0] == "A");
        REQUIRE(a2[0] == "G");

        REQUIRE(chrom[3] == "X");
        REQUIRE(pos[3] == 4000);
        REQUIRE(a1[3] == "T");
        REQUIRE(a2[3] == "C");
    }
}

TEST_CASE("read_bim - Invalid position data", "[reader][bim]")
{
    FileFixture files;

    SECTION("Exception - non-numeric position")
    {
        auto file_path
            = files.create_text_file("1\trs12345\t0\tinvalid\tA\tG", ".bim");
        REQUIRE_THROWS_AS(gelex::read_bim(file_path), gelex::GelexException);
    }
}

TEST_CASE("read_bim - Index lookup", "[reader][bim]")
{
    FileFixture files;

    auto file_path = files.create_text_file(
        "1\trs12345\t0\t1000\tA\tG\n"
        "2\trs67890\t0.001\t2000\tC\tT",
        ".bim");

    auto df = gelex::read_bim(file_path);

    REQUIRE(df.index().at("rs12345") == 0);
    REQUIRE(df.index().at("rs67890") == 1);
    REQUIRE(df.index().contains("rs12345"));
    REQUIRE_FALSE(df.index().contains("rs99999"));
}

TEST_CASE("read_bim - take_keys from rvalue", "[reader][bim]")
{
    FileFixture files;

    auto file_path = files.create_text_file(
        "1\trs12345\t0\t1000\tA\tG\n"
        "2\trs67890\t0.001\t2000\tC\tT",
        ".bim");

    auto ids = gelex::read_bim(file_path).index().take_keys();
    REQUIRE(ids.size() == 2);
    REQUIRE(ids[0] == "rs12345");
    REQUIRE(ids[1] == "rs67890");
}
