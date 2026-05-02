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
#include <span>
#include <string>
#include <vector>

#include <Eigen/Core>
#include <catch2/catch_test_macros.hpp>

#include "gelex/data/dataframe/column.h"
#include "gelex/exception.h"

using gelex::GelexException;
using gelex::dataframe::Column;

// ================================================================
// Construction
// ================================================================

TEST_CASE("Column name-only construction yields empty column", "[dataframe]")
{
    Column col("score");

    REQUIRE(col.name() == "score");
    REQUIRE(col.size() == 0);
}

TEST_CASE(
    "Column construction with int32_t vector stores correct data",
    "[dataframe]")
{
    std::vector<std::int32_t> data = {1, 2, 3};
    Column col("ids", data);

    REQUIRE(col.name() == "ids");
    REQUIRE(col.size() == 3);

    auto span = col.as<std::int32_t>();
    REQUIRE(span[0] == 1);
    REQUIRE(span[1] == 2);
    REQUIRE(span[2] == 3);
}

TEST_CASE(
    "Column construction with double vector stores correct data",
    "[dataframe]")
{
    std::vector<double> data = {1.1, 2.2, 3.3};
    Column col("vals", data);

    REQUIRE(col.name() == "vals");
    REQUIRE(col.size() == 3);

    auto span = col.as<double>();
    REQUIRE(span[0] == 1.1);
    REQUIRE(span[1] == 2.2);
    REQUIRE(span[2] == 3.3);
}

TEST_CASE(
    "Column construction with string vector stores correct data",
    "[dataframe]")
{
    std::vector<std::string> data = {"alpha", "beta", "gamma"};
    Column col("labels", data);

    REQUIRE(col.name() == "labels");
    REQUIRE(col.size() == 3);

    auto span = col.as<std::string>();
    REQUIRE(span[0] == "alpha");
    REQUIRE(span[1] == "beta");
    REQUIRE(span[2] == "gamma");
}

// ================================================================
// push_back
// ================================================================

TEST_CASE(
    "push_back int32_t into empty column initialises storage correctly",
    "[dataframe]")
{
    Column col("x");

    col.push_back(std::int32_t{42});

    REQUIRE(col.size() == 1);
    REQUIRE(col.as<std::int32_t>()[0] == 42);
}

TEST_CASE(
    "push_back double three times yields size 3 with correct values",
    "[dataframe]")
{
    Column col("d");

    col.push_back(1.0);
    col.push_back(2.5);
    col.push_back(3.14);

    REQUIRE(col.size() == 3);

    auto span = col.as<double>();
    REQUIRE(span[0] == 1.0);
    REQUIRE(span[1] == 2.5);
    REQUIRE(span[2] == 3.14);
}

TEST_CASE(
    "push_back string into empty column initialises storage correctly",
    "[dataframe]")
{
    Column col("s");

    col.push_back(std::string{"hello"});

    REQUIRE(col.size() == 1);
    REQUIRE(col.as<std::string>()[0] == "hello");
}

TEST_CASE(
    "push_back type mismatch after int init throws GelexException",
    "[dataframe]")
{
    Column col("mixed");

    col.push_back(std::int32_t{10});
    REQUIRE_THROWS_AS(col.push_back(std::string{"oops"}), GelexException);
}

// ================================================================
// as<T>()
// ================================================================

TEST_CASE("as<T> mutable overload returns modifiable span", "[dataframe]")
{
    Column col("v", std::vector<std::int32_t>{1, 2, 3});

    auto span = col.as<std::int32_t>();
    span[0] = 99;

    REQUIRE(col.as<std::int32_t>()[0] == 99);
}

TEST_CASE("as<T> const overload returns read-only span", "[dataframe]")
{
    const Column col("v", std::vector<double>{1.0, 2.0});

    auto span = col.as<double>();

    REQUIRE(span[0] == 1.0);
    REQUIRE(span[1] == 2.0);
}

TEST_CASE(
    "as<T> with wrong type on int column throws GelexException",
    "[dataframe]")
{
    Column col("c", std::vector<std::int32_t>{1, 2});

    REQUIRE_THROWS_AS(col.as<double>(), GelexException);
}

TEST_CASE(
    "as<T> const overload with wrong type throws GelexException",
    "[dataframe]")
{
    const Column col("c", std::vector<std::int32_t>{1, 2});

    REQUIRE_THROWS_AS(col.as<double>(), GelexException);
}

TEST_CASE(
    "as<T> on empty monostate column throws GelexException",
    "[dataframe]")
{
    Column col("empty");

    REQUIRE_THROWS_AS(col.as<std::int32_t>(), GelexException);
    REQUIRE_THROWS_AS(col.as<double>(), GelexException);
    REQUIRE_THROWS_AS(col.as<std::string>(), GelexException);
}

TEST_CASE(
    "as<T> const on empty monostate column throws GelexException",
    "[dataframe]")
{
    const Column col("empty");

    REQUIRE_THROWS_AS(col.as<std::int32_t>(), GelexException);
}

// ================================================================
// take
// ================================================================

TEST_CASE("take<int32_t> moves out underlying vector", "[dataframe]")
{
    Column col("ids", std::vector<std::int32_t>{1, 2, 3});

    auto taken = std::move(col).take<std::int32_t>();

    REQUIRE(taken.size() == 3);
    REQUIRE(taken[0] == 1);
    REQUIRE(taken[1] == 2);
    REQUIRE(taken[2] == 3);
}

TEST_CASE("take<double> moves out underlying vector", "[dataframe]")
{
    Column col("vals", std::vector<double>{1.1, 2.2, 3.3});

    auto taken = std::move(col).take<double>();

    REQUIRE(taken.size() == 3);
    REQUIRE(taken[0] == 1.1);
    REQUIRE(taken[1] == 2.2);
    REQUIRE(taken[2] == 3.3);
}

TEST_CASE("take<string> moves out underlying vector", "[dataframe]")
{
    Column col("labels", std::vector<std::string>{"alpha", "beta"});

    auto taken = std::move(col).take<std::string>();

    REQUIRE(taken.size() == 2);
    REQUIRE(taken[0] == "alpha");
    REQUIRE(taken[1] == "beta");
}

TEST_CASE("take with wrong type throws GelexException", "[dataframe]")
{
    Column col("c", std::vector<std::int32_t>{1, 2});

    REQUIRE_THROWS_AS(std::move(col).take<double>(), GelexException);
}

TEST_CASE("take on empty monostate column throws GelexException", "[dataframe]")
{
    Column col("empty");

    REQUIRE_THROWS_AS(std::move(col).take<std::int32_t>(), GelexException);
}

// ================================================================
// to_map
// ================================================================

TEST_CASE("to_map<double> returns view over column data", "[dataframe]")
{
    Column col("v", std::vector<double>{1.0, 2.0, 3.0});

    auto map = col.to_map<double>();

    REQUIRE(map.size() == 3);
    REQUIRE(map.data() == col.as<double>().data());
    REQUIRE(map.isApprox(Eigen::Vector3d{1.0, 2.0, 3.0}));
}

TEST_CASE("to_map<int32_t> returns view over int column data", "[dataframe]")
{
    Column col("i", std::vector<std::int32_t>{10, 20, 30});

    auto map = col.to_map<std::int32_t>();

    REQUIRE(map.size() == 3);
    REQUIRE(map.data() == col.as<std::int32_t>().data());
}

TEST_CASE("to_map with wrong type throws GelexException", "[dataframe]")
{
    Column col("v", std::vector<double>{1.0, 2.0});

    REQUIRE_THROWS_AS(col.to_map<std::int32_t>(), GelexException);
}

// ================================================================
// to_mat
// ================================================================

TEST_CASE("to_mat<double> returns copy of column data", "[dataframe]")
{
    Column col("v", std::vector<double>{1.0, 2.0, 3.0});

    auto vec = col.to_mat<double>();

    REQUIRE(vec.size() == 3);
    REQUIRE(vec.data() != col.as<double>().data());
    REQUIRE(vec.isApprox(Eigen::Vector3d{1.0, 2.0, 3.0}));
}

TEST_CASE("to_mat with wrong type throws GelexException", "[dataframe]")
{
    Column col("v", std::vector<double>{1.0});

    REQUIRE_THROWS_AS(col.to_mat<std::int32_t>(), GelexException);
}
