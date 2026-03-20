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

#include <catch2/catch_test_macros.hpp>

#include "gelex/data/dataframe/index.h"
#include "gelex/exception.h"

using gelex::DuplicateIndexException;
using gelex::InvalidInputException;
using gelex::df::Index;

// ----------------------------------------------------------------
// Helper: build an Index<string> from a brace-list
// ----------------------------------------------------------------
namespace
{

auto make_str_index(std::vector<std::string> keys) -> Index<std::string>
{
    return Index<std::string>(std::move(keys));
}

auto make_int_index(std::vector<std::int32_t> keys) -> Index<std::int32_t>
{
    return Index<std::int32_t>(std::move(keys));
}

}  // namespace

// ================================================================
// Construction
// ================================================================

TEST_CASE("Index<string> construction from vector", "[dataframe]")
{
    auto idx = make_str_index({"a", "b", "c"});

    REQUIRE(idx.size() == 3);
    REQUIRE(idx.contains("a"));
    REQUIRE(idx.contains("b"));
    REQUIRE(idx.contains("c"));
    REQUIRE(idx["a"] == 0);
    REQUIRE(idx["b"] == 1);
    REQUIRE(idx["c"] == 2);
}

TEST_CASE("Index<int32_t> construction from vector", "[dataframe]")
{
    auto idx = make_int_index({10, 20, 30});

    REQUIRE(idx.size() == 3);
    REQUIRE(idx.contains(10));
    REQUIRE(idx.contains(20));
    REQUIRE(idx.contains(30));
    REQUIRE(idx[10] == 0);
    REQUIRE(idx[20] == 1);
    REQUIRE(idx[30] == 2);
}

TEST_CASE(
    "Index construction from empty vector yields size zero",
    "[dataframe]")
{
    REQUIRE(Index<std::string>(std::vector<std::string>{}).size() == 0);
    REQUIRE(Index<std::int32_t>(std::vector<std::int32_t>{}).size() == 0);
}

TEST_CASE("Index<string> construction with duplicate key throws", "[dataframe]")
{
    REQUIRE_THROWS_AS(make_str_index({"a", "b", "a"}), DuplicateIndexException);
}

TEST_CASE(
    "Index<int32_t> construction with duplicate key throws",
    "[dataframe]")
{
    REQUIRE_THROWS_AS(make_int_index({10, 20, 10}), DuplicateIndexException);
}

// ================================================================
// push_back
// ================================================================

TEST_CASE(
    "Index<string> push_back increases size and assigns positions",
    "[dataframe]")
{
    Index<std::string> idx;

    idx.push_back("x");
    REQUIRE(idx.size() == 1);
    REQUIRE(idx["x"] == 0);

    idx.push_back("y");
    REQUIRE(idx.size() == 2);
    REQUIRE(idx["y"] == 1);

    idx.push_back("z");
    REQUIRE(idx.size() == 3);
    REQUIRE(idx["z"] == 2);
}

TEST_CASE(
    "Index<int32_t> push_back increases size and assigns positions",
    "[dataframe]")
{
    Index<std::int32_t> idx;

    idx.push_back(100);
    REQUIRE(idx.size() == 1);
    REQUIRE(idx[100] == 0);

    idx.push_back(200);
    REQUIRE(idx.size() == 2);
    REQUIRE(idx[200] == 1);

    idx.push_back(300);
    REQUIRE(idx.size() == 3);
    REQUIRE(idx[300] == 2);
}

TEST_CASE("Index<string> push_back of existing key throws", "[dataframe]")
{
    auto idx = make_str_index({"a", "b"});
    REQUIRE_THROWS_AS(idx.push_back("a"), DuplicateIndexException);
}

TEST_CASE("Index<int32_t> push_back of existing key throws", "[dataframe]")
{
    auto idx = make_int_index({10, 20});
    REQUIRE_THROWS_AS(idx.push_back(10), DuplicateIndexException);
}

// ================================================================
// operator[] and contains
// ================================================================

TEST_CASE("Index<string> operator[] returns correct position", "[dataframe]")
{
    auto idx = make_str_index({"foo", "bar", "baz"});
    REQUIRE(idx["foo"] == 0);
    REQUIRE(idx["bar"] == 1);
    REQUIRE(idx["baz"] == 2);
}

TEST_CASE("Index<string> operator[] for missing key throws", "[dataframe]")
{
    auto idx = make_str_index({"a", "b"});
    REQUIRE_THROWS_AS(idx["missing"], InvalidInputException);
}

TEST_CASE("Index<int32_t> operator[] for missing key throws", "[dataframe]")
{
    auto idx = make_int_index({1, 2, 3});
    REQUIRE_THROWS_AS(idx[99], InvalidInputException);
}

TEST_CASE(
    "Index<string> contains returns true and false correctly",
    "[dataframe]")
{
    auto idx = make_str_index({"p", "q"});
    REQUIRE(idx.contains("p"));
    REQUIRE(idx.contains("q"));
    REQUIRE_FALSE(idx.contains("r"));
}

TEST_CASE(
    "Index<int32_t> contains returns true and false correctly",
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

TEST_CASE("Index<string> data returns span over stored keys", "[dataframe]")
{
    auto idx = make_str_index({"a", "b", "c"});

    auto k = idx.data();

    REQUIRE(k.size() == 3);
    REQUIRE(k[0] == "a");
    REQUIRE(k[1] == "b");
    REQUIRE(k[2] == "c");
}

TEST_CASE("Index<int32_t> data returns span over stored keys", "[dataframe]")
{
    auto idx = make_int_index({10, 20, 30});

    auto k = idx.data();

    REQUIRE(k.size() == 3);
    REQUIRE(k[0] == 10);
    REQUIRE(k[1] == 20);
    REQUIRE(k[2] == 30);
}

TEST_CASE("Index<string> data on empty index returns empty span", "[dataframe]")
{
    Index<std::string> idx;
    REQUIRE(idx.data().empty());
}

// ================================================================
// take_keys
// ================================================================

TEST_CASE("Index<string> take_data moves out underlying vector", "[dataframe]")
{
    auto idx = make_str_index({"x", "y", "z"});

    auto taken = std::move(idx).take_data();

    REQUIRE(taken.size() == 3);
    REQUIRE(taken[0] == "x");
    REQUIRE(taken[1] == "y");
    REQUIRE(taken[2] == "z");
}

TEST_CASE("Index<int32_t> take_data moves out underlying vector", "[dataframe]")
{
    auto idx = make_int_index({100, 200, 300});

    auto taken = std::move(idx).take_data();

    REQUIRE(taken.size() == 3);
    REQUIRE(taken[0] == 100);
    REQUIRE(taken[1] == 200);
    REQUIRE(taken[2] == 300);
}

TEST_CASE(
    "Index<string> take_data on empty index returns empty vector",
    "[dataframe]")
{
    Index<std::string> idx;
    auto taken = std::move(idx).take_data();
    REQUIRE(taken.empty());
}

// ================================================================
// gather
// ================================================================

TEST_CASE(
    "Index<string> gather reorders keys and rebuilds positions",
    "[dataframe]")
{
    auto idx = make_str_index({"a", "b", "c", "d"});

    std::vector<std::size_t> indices = {2, 0, 3};
    idx.gather(indices);

    REQUIRE(idx.size() == 3);
    // After gather: [c, a, d]
    REQUIRE(idx["c"] == 0);
    REQUIRE(idx["a"] == 1);
    REQUIRE(idx["d"] == 2);
    REQUIRE_FALSE(idx.contains("b"));
}

TEST_CASE(
    "Index<int32_t> gather reorders keys and rebuilds positions",
    "[dataframe]")
{
    auto idx = make_int_index({10, 20, 30, 40});

    std::vector<std::size_t> indices = {3, 1};
    idx.gather(indices);

    REQUIRE(idx.size() == 2);
    REQUIRE(idx[40] == 0);
    REQUIRE(idx[20] == 1);
    REQUIRE_FALSE(idx.contains(10));
    REQUIRE_FALSE(idx.contains(30));
}

TEST_CASE(
    "Index<string> gather with empty indices yields size zero",
    "[dataframe]")
{
    auto idx = make_str_index({"a", "b", "c"});
    std::vector<std::size_t> empty_indices;
    idx.gather(empty_indices);
    REQUIRE(idx.size() == 0);
}

// ================================================================
// intersect
// ================================================================

TEST_CASE(
    "Index<string> intersect of two overlapping sets yields correct positions",
    "[dataframe]")
{
    auto idx0 = make_str_index({"a", "b", "c"});
    auto idx1 = make_str_index({"b", "c", "d"});

    // common keys sorted: {b, c}
    auto positions = Index<std::string>::intersect({&idx0, &idx1});

    REQUIRE(positions.size() == 2);
    // positions[0]: b→1, c→2 in idx0
    REQUIRE(positions[0] == std::vector<std::size_t>{1, 2});
    // positions[1]: b→0, c→1 in idx1
    REQUIRE(positions[1] == std::vector<std::size_t>{0, 1});
}

TEST_CASE(
    "Index<string> intersect of three overlapping sets yields correct "
    "positions",
    "[dataframe]")
{
    auto idx0 = make_str_index({"a", "b", "c"});
    auto idx1 = make_str_index({"b", "c", "d"});
    auto idx2 = make_str_index({"c", "d", "e"});

    // common key: {c} only
    auto positions = Index<std::string>::intersect({&idx0, &idx1, &idx2});

    REQUIRE(positions.size() == 3);
    REQUIRE(positions[0] == std::vector<std::size_t>{2});  // c→2 in idx0
    REQUIRE(positions[1] == std::vector<std::size_t>{1});  // c→1 in idx1
    REQUIRE(positions[2] == std::vector<std::size_t>{0});  // c→0 in idx2
}

TEST_CASE(
    "Index<string> intersect of disjoint sets returns empty position vectors",
    "[dataframe]")
{
    auto idx0 = make_str_index({"a", "b"});
    auto idx1 = make_str_index({"c", "d"});

    auto positions = Index<std::string>::intersect({&idx0, &idx1});

    REQUIRE(positions.size() == 2);
    REQUIRE(positions[0].empty());
    REQUIRE(positions[1].empty());
}

TEST_CASE(
    "Index<string> intersect with empty span returns empty result",
    "[dataframe]")
{
    std::span<const Index<std::string>*> empty_span;
    auto positions = Index<std::string>::intersect(empty_span);
    REQUIRE(positions.empty());
}

TEST_CASE(
    "Index<int32_t> intersect of two overlapping sets yields correct positions",
    "[dataframe]")
{
    auto idx0 = make_int_index({10, 20, 30});
    auto idx1 = make_int_index({20, 30, 40});

    // common keys sorted numerically: {20, 30}
    auto positions = Index<std::int32_t>::intersect({&idx0, &idx1});

    REQUIRE(positions.size() == 2);
    // positions[0]: 20→1, 30→2 in idx0
    REQUIRE(positions[0] == std::vector<std::size_t>{1, 2});
    // positions[1]: 20→0, 30→1 in idx1
    REQUIRE(positions[1] == std::vector<std::size_t>{0, 1});
}
