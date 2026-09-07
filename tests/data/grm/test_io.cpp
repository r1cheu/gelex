// Copyright 2026 RuLei Chen
// SPDX-License-Identifier: Apache-2.0

#include <Eigen/Core>
#include <catch2/catch_test_macros.hpp>
#include <string>
#include <vector>

#include "gelex/data/dataframe/index.h"
#include "gelex/data/grm/io.h"

#include "file_fixture.h"
#include "sample_id_fixture.h"

using gelex::test::FileFixture;

TEST_CASE("read_grm_ids reads GRM sample IDs", "[data][grm][io]")
{
    FileFixture files;
    auto prefix = files.generate_random_file_path("");
    std::vector<std::string> ids{
        gelex::make_sample_id("F1", "I1"),
        gelex::make_sample_id("F1", "I2"),
        gelex::make_sample_id("F2", "I3")};

    gelex::write_grm_ids(prefix.string(), ids);

    auto index = gelex::read_grm_ids(prefix.string());

    REQUIRE(index.size() == 3);
    REQUIRE(index.at(ids[0]) == 0);
    REQUIRE(index.at(ids[1]) == 1);
    REQUIRE(index.at(ids[2]) == 2);
}

TEST_CASE("read_grm reads full GRM matrix", "[data][grm][io]")
{
    FileFixture files;
    auto prefix = files.generate_random_file_path("");
    std::vector<std::string> ids{
        gelex::make_sample_id("F1", "I1"),
        gelex::make_sample_id("F1", "I2"),
        gelex::make_sample_id("F2", "I3")};
    Eigen::MatrixXd matrix{
        {1.0, 0.25, 0.5}, {0.25, 2.0, 0.75}, {0.5, 0.75, 4.0}};

    gelex::write_grm(prefix.string(), matrix, ids);

    auto actual = gelex::read_grm(prefix.string(), nullptr, false);

    REQUIRE(actual.isApprox(matrix));
}

TEST_CASE("read_grm normalizes by mean diagonal", "[data][grm][io]")
{
    FileFixture files;
    auto prefix = files.generate_random_file_path("");
    std::vector<std::string> ids{
        gelex::make_sample_id("F1", "I1"),
        gelex::make_sample_id("F1", "I2"),
        gelex::make_sample_id("F2", "I3")};
    Eigen::MatrixXd matrix{
        {1.0, 0.25, 0.5}, {0.25, 2.0, 0.75}, {0.5, 0.75, 3.0}};

    gelex::write_grm(prefix.string(), matrix, ids);

    auto actual = gelex::read_grm(prefix.string());
    Eigen::MatrixXd expected = matrix / 2.0;

    REQUIRE(actual.isApprox(expected));
}

TEST_CASE("read_grm reads reordered GRM subset", "[data][grm][io]")
{
    FileFixture files;
    auto prefix = files.generate_random_file_path("");
    std::vector<std::string> ids{
        gelex::make_sample_id("F1", "I1"),
        gelex::make_sample_id("F1", "I2"),
        gelex::make_sample_id("F2", "I3")};
    Eigen::MatrixXd matrix{
        {1.0, 0.25, 0.5}, {0.25, 2.0, 0.75}, {0.5, 0.75, 3.0}};

    gelex::write_grm(prefix.string(), matrix, ids);

    gelex::DataFrameIndex<std::string> target_index(
        std::vector<std::string>{ids[2], ids[0]});

    auto actual = gelex::read_grm(prefix.string(), &target_index, false);
    Eigen::MatrixXd expected{{3.0, 0.5}, {0.5, 1.0}};

    REQUIRE(actual.isApprox(expected));
}
