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
#include <utility>
#include <vector>

#include <catch2/catch_test_macros.hpp>

#include "gelex/data/dataframe/index.h"
#include "gelex/exception.h"

using gelex::DataFrameIndex;
using gelex::GelexException;
using gelex::intersect;

// ----------------------------------------------------------------
// Helper: build an DataFrameIndex<string> from a brace-list
// ----------------------------------------------------------------
namespace
{

auto make_str_index(std::vector<std::string> keys)
    -> DataFrameIndex<std::string>
{
    return DataFrameIndex<std::string>(std::move(keys));
}

auto make_int_index(std::vector<std::int32_t> keys)
    -> DataFrameIndex<std::int32_t>
{
    return DataFrameIndex<std::int32_t>(std::move(keys));
}

}  // namespace

// ================================================================
// Construction
// ================================================================

TEST_CASE("DataFrameIndex<string> construction from vector", "[dataframe]")
{
    auto idx = make_str_index({"a", "b", "c"});

    REQUIRE(idx.size() == 3);
    REQUIRE(idx.contains("a"));
    REQUIRE(idx.contains("b"));
    REQUIRE(idx.contains("c"));
    REQUIRE(idx.at("a") == 0);
    REQUIRE(idx.at("b") == 1);
    REQUIRE(idx.at("c") == 2);
}

TEST_CASE("DataFrameIndex<int32_t> construction from vector", "[dataframe]")
{
    auto idx = make_int_index({10, 20, 30});

    REQUIRE(idx.size() == 3);
    REQUIRE(idx.contains(10));
    REQUIRE(idx.contains(20));
    REQUIRE(idx.contains(30));
    REQUIRE(idx.at(10) == 0);
    REQUIRE(idx.at(20) == 1);
    REQUIRE(idx.at(30) == 2);
}

TEST_CASE(
    "DataFrameIndex construction from empty vector yields size zero",
    "[dataframe]")
{
    REQUIRE(
        DataFrameIndex<std::string>(std::vector<std::string>{}).size() == 0);
    REQUIRE(
        DataFrameIndex<std::int32_t>(std::vector<std::int32_t>{}).size() == 0);
}

TEST_CASE(
    "DataFrameIndex<string> construction with duplicate key throws",
    "[dataframe]")
{
    REQUIRE_THROWS_AS(make_str_index({"a", "b", "a"}), GelexException);
}

TEST_CASE(
    "DataFrameIndex<int32_t> construction with duplicate key throws",
    "[dataframe]")
{
    REQUIRE_THROWS_AS(make_int_index({10, 20, 10}), GelexException);
}

// ================================================================
// push_back
// ================================================================

TEST_CASE(
    "DataFrameIndex<string> push_back increases size and assigns positions",
    "[dataframe]")
{
    DataFrameIndex<std::string> idx;

    idx.push_back("x");
    REQUIRE(idx.size() == 1);
    REQUIRE(idx.at("x") == 0);

    idx.push_back("y");
    REQUIRE(idx.size() == 2);
    REQUIRE(idx.at("y") == 1);

    idx.push_back("z");
    REQUIRE(idx.size() == 3);
    REQUIRE(idx.at("z") == 2);
}

TEST_CASE(
    "DataFrameIndex<int32_t> push_back increases size and assigns positions",
    "[dataframe]")
{
    DataFrameIndex<std::int32_t> idx;

    idx.push_back(100);
    REQUIRE(idx.size() == 1);
    REQUIRE(idx.at(100) == 0);

    idx.push_back(200);
    REQUIRE(idx.size() == 2);
    REQUIRE(idx.at(200) == 1);

    idx.push_back(300);
    REQUIRE(idx.size() == 3);
    REQUIRE(idx.at(300) == 2);
}

TEST_CASE(
    "DataFrameIndex<string> push_back of existing key throws",
    "[dataframe]")
{
    auto idx = make_str_index({"a", "b"});
    REQUIRE_THROWS_AS(idx.push_back("a"), GelexException);
}

TEST_CASE(
    "DataFrameIndex<int32_t> push_back of existing key throws",
    "[dataframe]")
{
    auto idx = make_int_index({10, 20});
    REQUIRE_THROWS_AS(idx.push_back(10), GelexException);
}

// ================================================================
// operator[] and contains
// ================================================================

TEST_CASE(
    "DataFrameIndex<string> operator[] returns correct position",
    "[dataframe]")
{
    auto idx = make_str_index({"foo", "bar", "baz"});
    REQUIRE(idx.at("foo") == 0);
    REQUIRE(idx.at("bar") == 1);
    REQUIRE(idx.at("baz") == 2);
}

TEST_CASE(
    "DataFrameIndex<string> operator[] for missing key throws",
    "[dataframe]")
{
    auto idx = make_str_index({"a", "b"});
    REQUIRE_THROWS_AS(idx.at("missing"), GelexException);
}

TEST_CASE(
    "DataFrameIndex<int32_t> operator[] for missing key throws",
    "[dataframe]")
{
    auto idx = make_int_index({1, 2, 3});
    REQUIRE_THROWS_AS(idx.at(99), GelexException);
}

TEST_CASE(
    "DataFrameIndex<string> contains returns true and false correctly",
    "[dataframe]")
{
    auto idx = make_str_index({"p", "q"});
    REQUIRE(idx.contains("p"));
    REQUIRE(idx.contains("q"));
    REQUIRE_FALSE(idx.contains("r"));
}

TEST_CASE(
    "DataFrameIndex<int32_t> contains returns true and false correctly",
    "[dataframe]")
{
    auto idx = make_int_index({7, 42});
    REQUIRE(idx.contains(7));
    REQUIRE(idx.contains(42));
    REQUIRE_FALSE(idx.contains(0));
}

// ================================================================
// keys
// ================================================================

TEST_CASE(
    "DataFrameIndex<string> data returns span over stored keys",
    "[dataframe]")
{
    auto idx = make_str_index({"a", "b", "c"});

    auto k = idx.keys();

    REQUIRE(k.size() == 3);
    REQUIRE(k[0] == "a");
    REQUIRE(k[1] == "b");
    REQUIRE(k[2] == "c");
}

TEST_CASE(
    "DataFrameIndex<int32_t> data returns span over stored keys",
    "[dataframe]")
{
    auto idx = make_int_index({10, 20, 30});

    auto k = idx.keys();

    REQUIRE(k.size() == 3);
    REQUIRE(k[0] == 10);
    REQUIRE(k[1] == 20);
    REQUIRE(k[2] == 30);
}

TEST_CASE(
    "DataFrameIndex<string> data on empty index returns empty span",
    "[dataframe]")
{
    DataFrameIndex<std::string> idx;
    REQUIRE(idx.keys().empty());
}

// ================================================================
// take_keys
// ================================================================

TEST_CASE(
    "DataFrameIndex<string> take_data moves out underlying vector",
    "[dataframe]")
{
    auto idx = make_str_index({"x", "y", "z"});

    auto taken = std::move(idx).take_keys();

    REQUIRE(taken.size() == 3);
    REQUIRE(taken[0] == "x");
    REQUIRE(taken[1] == "y");
    REQUIRE(taken[2] == "z");
}

TEST_CASE(
    "DataFrameIndex<int32_t> take_data moves out underlying vector",
    "[dataframe]")
{
    auto idx = make_int_index({100, 200, 300});

    auto taken = std::move(idx).take_keys();

    REQUIRE(taken.size() == 3);
    REQUIRE(taken[0] == 100);
    REQUIRE(taken[1] == 200);
    REQUIRE(taken[2] == 300);
}

TEST_CASE(
    "DataFrameIndex<string> take_data on empty index returns empty vector",
    "[dataframe]")
{
    DataFrameIndex<std::string> idx;
    auto taken = std::move(idx).take_keys();
    REQUIRE(taken.empty());
}

// ================================================================
// gather
// ================================================================

TEST_CASE(
    "DataFrameIndex<string> gather reorders keys and rebuilds positions",
    "[dataframe]")
{
    auto idx = make_str_index({"a", "b", "c", "d"});

    // target: [c, a, d]
    idx.gather(make_str_index({"c", "a", "d"}));

    REQUIRE(idx.size() == 3);
    REQUIRE(idx.at("c") == 0);
    REQUIRE(idx.at("a") == 1);
    REQUIRE(idx.at("d") == 2);
    REQUIRE_FALSE(idx.contains("b"));
}

TEST_CASE(
    "DataFrameIndex<int32_t> gather reorders keys and rebuilds positions",
    "[dataframe]")
{
    auto idx = make_int_index({10, 20, 30, 40});

    // target: [40, 20]
    idx.gather(make_int_index({40, 20}));

    REQUIRE(idx.size() == 2);
    REQUIRE(idx.at(40) == 0);
    REQUIRE(idx.at(20) == 1);
    REQUIRE_FALSE(idx.contains(10));
    REQUIRE_FALSE(idx.contains(30));
}

TEST_CASE(
    "DataFrameIndex<string> gather with empty index yields size zero",
    "[dataframe]")
{
    auto idx = make_str_index({"a", "b", "c"});
    idx.gather(DataFrameIndex<std::string>{});
    REQUIRE(idx.size() == 0);
}

// ================================================================
// intersect
// ================================================================

TEST_CASE(
    "DataFrameIndex<string> intersect of two overlapping sets returns common "
    "keys",
    "[dataframe]")
{
    auto idx0 = make_str_index({"a", "b", "c"});
    auto idx1 = make_str_index({"b", "c", "d"});

    auto common = intersect({&idx0, &idx1});

    REQUIRE(common.size() == 2);
    REQUIRE(common.keys()[0] == "b");
    REQUIRE(common.keys()[1] == "c");
}

TEST_CASE(
    "DataFrameIndex<string> intersect of three overlapping sets returns common "
    "keys",
    "[dataframe]")
{
    auto idx0 = make_str_index({"a", "b", "c"});
    auto idx1 = make_str_index({"b", "c", "d"});
    auto idx2 = make_str_index({"c", "d", "e"});

    auto common = intersect({&idx0, &idx1, &idx2});

    REQUIRE(common.size() == 1);
    REQUIRE(common.keys()[0] == "c");
}

TEST_CASE(
    "DataFrameIndex<string> intersect of disjoint sets returns empty "
    "DataFrameIndex",
    "[dataframe]")
{
    auto idx0 = make_str_index({"a", "b"});
    auto idx1 = make_str_index({"c", "d"});

    auto common = intersect({&idx0, &idx1});

    REQUIRE(common.size() == 0);
}

TEST_CASE(
    "DataFrameIndex<string> intersect with empty span returns empty "
    "DataFrameIndex",
    "[dataframe]")
{
    std::span<const DataFrameIndex<std::string>* const> empty_span;
    auto common = intersect(empty_span);
    REQUIRE(common.size() == 0);
}

TEST_CASE(
    "DataFrameIndex<int32_t> intersect of two overlapping sets returns common "
    "keys",
    "[dataframe]")
{
    auto idx0 = make_int_index({10, 20, 30});
    auto idx1 = make_int_index({20, 30, 40});

    auto common = intersect({&idx0, &idx1});

    REQUIRE(common.size() == 2);
    REQUIRE(common.keys()[0] == 20);
    REQUIRE(common.keys()[1] == 30);
}
