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
#include "gelex/data/frame/dataframe.h"
#include "gelex/data/frame/dummy_encode.h"

using gelex::DataFrame;
using gelex::DummyEncode;
using gelex::test::FileFixture;

TEST_CASE("DummyEncode builds matrix and metadata", "[data][dummy_encode]")
{
    FileFixture files;
    auto path = files.create_text_file(
        "FID\tIID\tGroup\tSite\n"
        "f1\ti1\tB\tX\n"
        "f2\ti2\tA\tY\n"
        "f3\ti3\tC\tX\n"
        "f4\ti4\tA\tZ\n");

    auto frame = DataFrame<std::string>::read(path);
    auto dcov = DummyEncode(frame);

    REQUIRE(dcov.names.size() == 2);
    REQUIRE(dcov.names[0] == "Group");
    REQUIRE(dcov.names[1] == "Site");

    REQUIRE(dcov.levels.size() == 2);
    REQUIRE(dcov.levels[0].size() == 3);
    REQUIRE(dcov.levels[0][0] == "B");
    REQUIRE(dcov.levels[0][1] == "A");
    REQUIRE(dcov.levels[0][2] == "C");

    REQUIRE(dcov.levels[1].size() == 3);
    REQUIRE(dcov.levels[1][0] == "X");
    REQUIRE(dcov.levels[1][1] == "Y");
    REQUIRE(dcov.levels[1][2] == "Z");

    REQUIRE(dcov.reference_levels.size() == 2);
    REQUIRE(dcov.reference_levels[0] == "B");
    REQUIRE(dcov.reference_levels[1] == "X");

    REQUIRE(dcov.X.rows() == 4);
    REQUIRE(dcov.X.cols() == 4);

    REQUIRE(dcov.X(0, 0) == 0.0);  // Group=A
    REQUIRE(dcov.X(0, 1) == 0.0);  // Group=C
    REQUIRE(dcov.X(0, 2) == 0.0);  // Site=Y
    REQUIRE(dcov.X(0, 3) == 0.0);  // Site=Z

    REQUIRE(dcov.X(1, 0) == 1.0);  // Group=A
    REQUIRE(dcov.X(1, 1) == 0.0);  // Group=C
    REQUIRE(dcov.X(1, 2) == 1.0);  // Site=Y
    REQUIRE(dcov.X(1, 3) == 0.0);  // Site=Z

    REQUIRE(dcov.X(2, 0) == 0.0);  // Group=A
    REQUIRE(dcov.X(2, 1) == 1.0);  // Group=C
    REQUIRE(dcov.X(2, 2) == 0.0);  // Site=Y
    REQUIRE(dcov.X(2, 3) == 0.0);  // Site=Z

    REQUIRE(dcov.X(3, 0) == 1.0);  // Group=A
    REQUIRE(dcov.X(3, 1) == 0.0);  // Group=C
    REQUIRE(dcov.X(3, 2) == 0.0);  // Site=Y
    REQUIRE(dcov.X(3, 3) == 1.0);  // Site=Z
}

TEST_CASE(
    "DummyEncode skips single-level columns",
    "[data][dummy_encode]")
{
    FileFixture files;
    auto path = files.create_text_file(
        "FID\tIID\tConstant\tGroup\n"
        "f1\ti1\tK\tA\n"
        "f2\ti2\tK\tB\n"
        "f3\ti3\tK\tA\n");

    auto frame = DataFrame<std::string>::read(path);
    auto dcov = DummyEncode(frame);

    REQUIRE(dcov.names.size() == 1);
    REQUIRE(dcov.names[0] == "Group");

    REQUIRE(dcov.levels.size() == 1);
    REQUIRE(dcov.reference_levels.size() == 1);
    REQUIRE(dcov.reference_levels[0] == "A");

    REQUIRE(dcov.X.rows() == 3);
    REQUIRE(dcov.X.cols() == 1);
    REQUIRE(dcov.X(0, 0) == 0.0);
    REQUIRE(dcov.X(1, 0) == 1.0);
    REQUIRE(dcov.X(2, 0) == 0.0);
}

TEST_CASE(
    "DummyEncode returns empty matrix when all columns single-level",
    "[data][dummy_encode]")
{
    FileFixture files;
    auto path = files.create_text_file(
        "FID\tIID\tGroup\tSite\n"
        "f1\ti1\tA\tX\n"
        "f2\ti2\tA\tX\n");

    auto frame = DataFrame<std::string>::read(path);
    auto dcov = DummyEncode(frame);

    REQUIRE(dcov.names.empty());
    REQUIRE(dcov.levels.empty());
    REQUIRE(dcov.reference_levels.empty());
    REQUIRE(dcov.X.rows() == 2);
    REQUIRE(dcov.X.cols() == 0);
}
