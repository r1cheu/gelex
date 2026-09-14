// Copyright 2026 RuLei Chen
// SPDX-License-Identifier: Apache-2.0

#include <catch2/catch_test_macros.hpp>
#include <fstream>
#include <string>
#include <vector>

#include "gelex/data/sample_id_io.h"
#include "gelex/exception.h"

#include "file_fixture.h"
#include "sample_id_fixture.h"

using gelex::test::FileFixture;

TEST_CASE("sample IDs round-trip as <prefix>.id", "[data][sample_id]")
{
    FileFixture files;
    const auto prefix = files.generate_random_file_path("").string();
    const std::vector<std::string> ids{
        gelex::make_sample_id("F1", "I1"),
        gelex::make_sample_id("F1", "I2"),
        gelex::make_sample_id("F2", "I3")};

    gelex::write_sample_ids(prefix, ids);

    std::ifstream in{prefix + ".id"};
    std::string line;
    std::getline(in, line);
    REQUIRE(line == "F1\tI1");

    const auto index = gelex::read_sample_ids(prefix);
    REQUIRE(index.size() == 3);
    REQUIRE(index.at(ids[0]) == 0);
    REQUIRE(index.at(ids[2]) == 2);
}

TEST_CASE("read_sample_ids rejects an empty file", "[data][sample_id]")
{
    FileFixture files;
    const auto prefix = files.generate_random_file_path("").string();
    gelex::write_sample_ids(prefix, {});
    REQUIRE_THROWS_AS(gelex::read_sample_ids(prefix), gelex::GelexException);
}
