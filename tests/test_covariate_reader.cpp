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
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include "file_fixture.h"
#include "gelex/exception.h"
#include "gelex/predict/covariate_reader.h"

using Catch::Matchers::WithinAbs;
using gelex::CovariateParams;
using gelex::CovariateReader;
using gelex::test::FileFixture;

TEST_CASE("CovariateReader - normal file", "[predict][covariate_reader]")
{
    FileFixture files;
    auto path = files.create_text_file(
        "term\tmean\tstddev\n"
        "Intercept\t2.5\t0.1\n"
        "Age\t0.5\t0.05\n"
        "Sex_M\t-0.3\t0.02\n");

    CovariateReader reader(path);
    const auto& p = reader.params();

    REQUIRE(p.term_names.size() == 3);
    REQUIRE(p.coefficients.size() == 3);
    REQUIRE(p.term_names[0] == "Intercept");
    REQUIRE(p.term_names[1] == "Age");
    REQUIRE(p.term_names[2] == "Sex_M");
    REQUIRE_THAT(p.coefficients(0), WithinAbs(2.5, 1e-12));
    REQUIRE_THAT(p.coefficients(1), WithinAbs(0.5, 1e-12));
    REQUIRE_THAT(p.coefficients(2), WithinAbs(-0.3, 1e-12));
}

TEST_CASE("CovariateReader - single term", "[predict][covariate_reader]")
{
    FileFixture files;
    auto path = files.create_text_file(
        "term\tmean\tstddev\n"
        "Intercept\t1.8\t0.2\n");

    CovariateReader reader(path);
    REQUIRE(reader.params().term_names.size() == 1);
    REQUIRE(reader.params().term_names[0] == "Intercept");
    REQUIRE_THAT(reader.params().coefficients(0), WithinAbs(1.8, 1e-12));
}

TEST_CASE(
    "CovariateReader - extra columns ignored",
    "[predict][covariate_reader]")
{
    FileFixture files;
    auto path = files.create_text_file(
        "term\tmean\tstddev\tpercentile_5\tpercentile_95\tess\trhat\n"
        "Intercept\t2.5\t0.1\t2.3\t2.7\t1000\t1.0\n"
        "Age\t0.5\t0.05\t0.4\t0.6\t800\t1.01\n");

    CovariateReader reader(path);
    REQUIRE(reader.params().term_names.size() == 2);
    REQUIRE_THAT(reader.params().coefficients(0), WithinAbs(2.5, 1e-12));
    REQUIRE_THAT(reader.params().coefficients(1), WithinAbs(0.5, 1e-12));
}

TEST_CASE("CovariateReader - empty file throws", "[predict][covariate_reader]")
{
    FileFixture files;
    auto path = files.create_empty_file(".txt");
    REQUIRE_THROWS_AS(CovariateReader(path), gelex::FileFormatException);
}

TEST_CASE("CovariateReader - header only throws", "[predict][covariate_reader]")
{
    FileFixture files;
    auto path = files.create_text_file("term\tmean\tstddev\n");
    REQUIRE_THROWS_AS(CovariateReader(path), gelex::FileFormatException);
}

TEST_CASE(
    "CovariateReader - malformed lines skipped",
    "[predict][covariate_reader]")
{
    FileFixture files;
    auto path = files.create_text_file(
        "term\tmean\tstddev\n"
        "Intercept\t2.0\t0.1\n"
        "InvalidLine\n"
        "Age\t0.5\t0.05\n");

    CovariateReader reader(path);
    const auto& p = reader.params();
    REQUIRE(p.term_names.size() == 2);
    REQUIRE(p.term_names[0] == "Intercept");
    REQUIRE(p.term_names[1] == "Age");
}

TEST_CASE(
    "CovariateReader - non-numeric mean skipped",
    "[predict][covariate_reader]")
{
    FileFixture files;
    auto path = files.create_text_file(
        "term\tmean\tstddev\n"
        "Intercept\t1.0\t0.1\n"
        "Age\tnot_a_number\t0.05\n"
        "Height\t0.2\t0.03\n");

    CovariateReader reader(path);
    const auto& p = reader.params();
    REQUIRE(p.term_names.size() == 2);
    REQUIRE(p.term_names[0] == "Intercept");
    REQUIRE(p.term_names[1] == "Height");
}

TEST_CASE(
    "CovariateReader - empty mean field skipped",
    "[predict][covariate_reader]")
{
    FileFixture files;
    auto path = files.create_text_file(
        "term\tmean\tstddev\n"
        "Intercept\t1.0\t0.1\n"
        "Age\t\t0.05\n"
        "Height\t0.2\t0.03\n");

    CovariateReader reader(path);
    REQUIRE(reader.params().term_names.size() == 2);
}

TEST_CASE("CovariateReader - file not found", "[predict][covariate_reader]")
{
    FileFixture files;
    auto path = files.generate_random_file_path(".txt");
    REQUIRE_THROWS(CovariateReader(path));
}

TEST_CASE(
    "CovariateReader - duplicate terms kept in order",
    "[predict][covariate_reader]")
{
    FileFixture files;
    auto path = files.create_text_file(
        "term\tmean\tstddev\n"
        "Intercept\t1.0\t0.1\n"
        "Age\t0.3\t0.05\n"
        "Age\t0.5\t0.05\n");

    CovariateReader reader(path);
    const auto& p = reader.params();
    REQUIRE(p.term_names.size() == 3);
    REQUIRE(p.term_names[1] == "Age");
    REQUIRE(p.term_names[2] == "Age");
    REQUIRE_THAT(p.coefficients(1), WithinAbs(0.3, 1e-12));
    REQUIRE_THAT(p.coefficients(2), WithinAbs(0.5, 1e-12));
}

TEST_CASE(
    "CovariateReader - params() returns const ref",
    "[predict][covariate_reader]")
{
    FileFixture files;
    auto path = files.create_text_file(
        "term\tmean\tstddev\n"
        "Intercept\t1.0\t0.1\n");

    CovariateReader reader(path);
    const auto& p1 = reader.params();
    const auto& p2 = reader.params();
    REQUIRE(&p1 == &p2);
}

TEST_CASE(
    "CovariateReader - take_params() moves data",
    "[predict][covariate_reader]")
{
    FileFixture files;
    auto path = files.create_text_file(
        "term\tmean\tstddev\n"
        "Intercept\t1.0\t0.1\n"
        "Age\t0.5\t0.05\n");

    CovariateReader reader(path);
    REQUIRE(reader.params().term_names.size() == 2);

    CovariateParams moved = std::move(reader).take_params();
    REQUIRE(moved.term_names.size() == 2);
    REQUIRE(reader.params().term_names.empty());
}
