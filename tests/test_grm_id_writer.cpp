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

#include <algorithm>
#include <filesystem>
#include <fstream>
#include <sstream>
#include <string>
#include <string_view>
#include <vector>

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_exception.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>

#include "../src/data/grm_id_writer.h"
#include "file_fixture.h"
#include "gelex/data/dataframe_policy.h"
#include "gelex/exception.h"

namespace fs = std::filesystem;

using namespace gelex::detail;  // NOLINT
using Catch::Matchers::ContainsSubstring;
using Catch::Matchers::MessageMatches;
using gelex::test::FileFixture;

namespace
{

auto sid(std::string_view fid, std::string_view iid) -> std::string
{
    return gelex::make_sample_id(fid, iid);
}

auto read_file_content(const fs::path& file_path) -> std::string
{
    std::ifstream file(file_path);
    std::stringstream buffer;
    buffer << file.rdbuf();
    return buffer.str();
}

auto read_file_lines(const fs::path& file_path) -> std::vector<std::string>
{
    std::ifstream file(file_path);
    std::vector<std::string> lines;
    std::string line;
    while (std::getline(file, line))
    {
        lines.push_back(line);
    }
    return lines;
}

}  // namespace

TEST_CASE(
    "GrmIdWriter - Constructor and path access",
    "[grm_id_writer][construction]")
{
    FileFixture files;

    SECTION("Happy path - create writer with valid path")
    {
        auto file_path = files.generate_random_file_path(".grm.id");

        REQUIRE_NOTHROW(
            [&]()
            {
                GrmIdWriter writer(file_path);
                REQUIRE(writer.path() == file_path);
            }());
    }

    SECTION("Happy path - file is created on construction")
    {
        auto file_path = files.generate_random_file_path(".grm.id");

        {
            GrmIdWriter writer(file_path);
        }

        REQUIRE(fs::exists(file_path));
    }

    SECTION("Happy path - constructor truncates existing file")
    {
        auto file_path = files.generate_random_file_path(".grm.id");
        {
            std::ofstream out(file_path);
            out << "legacy-content\n";
        }
        REQUIRE(fs::file_size(file_path) > 0);

        {
            GrmIdWriter writer(file_path);
        }

        REQUIRE(fs::exists(file_path));
        REQUIRE(fs::file_size(file_path) == 0);
    }
}

TEST_CASE("GrmIdWriter - Write empty ID list", "[grm_id_writer][empty]")
{
    FileFixture files;
    auto file_path = files.generate_random_file_path(".grm.id");
    std::vector<std::string> ids;

    REQUIRE_NOTHROW(
        [&]()
        {
            GrmIdWriter writer(file_path);
            writer.write(ids);
        }());

    REQUIRE(fs::exists(file_path));
    REQUIRE(fs::file_size(file_path) == 0);
}

TEST_CASE("GrmIdWriter - Write canonical IDs", "[grm_id_writer][basic]")
{
    FileFixture files;

    SECTION("Happy path - write single canonical ID")
    {
        auto file_path = files.generate_random_file_path(".grm.id");
        std::vector<std::string> ids = {sid("FAM1", "IND1")};

        {
            GrmIdWriter writer(file_path);
            writer.write(ids);
        }

        REQUIRE(read_file_content(file_path) == "FAM1\tIND1\n");
    }

    SECTION("Happy path - write multiple canonical IDs")
    {
        auto file_path = files.generate_random_file_path(".grm.id");
        std::vector<std::string> ids
            = {sid("FAM1", "IND1"), sid("FAM2", "IND2"), sid("FAM3", "IND3")};

        {
            GrmIdWriter writer(file_path);
            writer.write(ids);
        }

        auto lines = read_file_lines(file_path);
        REQUIRE(lines.size() == 3);
        REQUIRE(lines[0] == "FAM1\tIND1");
        REQUIRE(lines[1] == "FAM2\tIND2");
        REQUIRE(lines[2] == "FAM3\tIND3");
    }

    SECTION("Happy path - FID and IID may include underscores")
    {
        auto file_path = files.generate_random_file_path(".grm.id");
        std::vector<std::string> ids
            = {sid("FAM_1", "IND_2"), sid("A_B_C", "X_Y_Z")};

        {
            GrmIdWriter writer(file_path);
            writer.write(ids);
        }

        auto lines = read_file_lines(file_path);
        REQUIRE(lines.size() == 2);
        REQUIRE(lines[0] == "FAM_1\tIND_2");
        REQUIRE(lines[1] == "A_B_C\tX_Y_Z");
    }
}

TEST_CASE("GrmIdWriter - Reject non-canonical IDs", "[grm_id_writer][error]")
{
    FileFixture files;

    SECTION("Exception - reject legacy FID_IID ID")
    {
        auto file_path = files.generate_random_file_path(".grm.id");
        std::vector<std::string> ids = {"FAM1_IND1"};

        GrmIdWriter writer(file_path);
        REQUIRE_THROWS_MATCHES(
            writer.write(ids),
            gelex::ArgumentValidationException,
            MessageMatches(ContainsSubstring("canonical FID<US>IID format")));
    }

    SECTION("Exception - reject ID without separator")
    {
        auto file_path = files.generate_random_file_path(".grm.id");
        std::vector<std::string> ids = {"NOSPLIT"};

        GrmIdWriter writer(file_path);
        REQUIRE_THROWS_MATCHES(
            writer.write(ids),
            gelex::ArgumentValidationException,
            MessageMatches(ContainsSubstring("canonical FID<US>IID format")));
    }

    SECTION("Exception - reject empty ID")
    {
        auto file_path = files.generate_random_file_path(".grm.id");
        std::vector<std::string> ids = {""};

        GrmIdWriter writer(file_path);
        REQUIRE_THROWS_MATCHES(
            writer.write(ids),
            gelex::ArgumentValidationException,
            MessageMatches(ContainsSubstring("cannot be empty")));
    }

    SECTION("Exception - reject empty FID")
    {
        auto file_path = files.generate_random_file_path(".grm.id");
        std::vector<std::string> ids
            = {std::string(1, gelex::kSampleIdSeparator) + "IND1"};

        GrmIdWriter writer(file_path);
        REQUIRE_THROWS_MATCHES(
            writer.write(ids),
            gelex::ArgumentValidationException,
            MessageMatches(ContainsSubstring("non-empty FID and IID")));
    }

    SECTION("Exception - reject empty IID")
    {
        auto file_path = files.generate_random_file_path(".grm.id");
        std::vector<std::string> ids
            = {"FAM1" + std::string(1, gelex::kSampleIdSeparator)};

        GrmIdWriter writer(file_path);
        REQUIRE_THROWS_MATCHES(
            writer.write(ids),
            gelex::ArgumentValidationException,
            MessageMatches(ContainsSubstring("non-empty FID and IID")));
    }
}

TEST_CASE("GrmIdWriter - Output format verification", "[grm_id_writer][format]")
{
    FileFixture files;

    SECTION("Happy path - verify tab separator and newline")
    {
        auto file_path = files.generate_random_file_path(".grm.id");
        std::vector<std::string> ids = {sid("A", "B"), sid("C", "D")};

        {
            GrmIdWriter writer(file_path);
            writer.write(ids);
        }

        auto content = read_file_content(file_path);
        REQUIRE(content == "A\tB\nC\tD\n");

        auto tab_count = std::count(content.begin(), content.end(), '\t');
        auto newline_count = std::count(content.begin(), content.end(), '\n');
        REQUIRE(tab_count == 2);
        REQUIRE(newline_count == 2);
    }

    SECTION("Happy path - file ends with newline")
    {
        auto file_path = files.generate_random_file_path(".grm.id");
        std::vector<std::string> ids = {sid("X", "Y")};

        {
            GrmIdWriter writer(file_path);
            writer.write(ids);
        }

        auto content = read_file_content(file_path);
        REQUIRE(!content.empty());
        REQUIRE(content.back() == '\n');
    }
}

TEST_CASE("GrmIdWriter - Multiple write calls", "[grm_id_writer][multiple]")
{
    FileFixture files;

    SECTION("Happy path - multiple write calls append")
    {
        auto file_path = files.generate_random_file_path(".grm.id");
        std::vector<std::string> ids1
            = {sid("FAM1", "IND1"), sid("FAM2", "IND2")};
        std::vector<std::string> ids2 = {sid("FAM3", "IND3")};
        std::vector<std::string> ids3
            = {sid("FAM4", "IND4"), sid("FAM5", "IND5")};

        {
            GrmIdWriter writer(file_path);
            writer.write(ids1);
            writer.write(ids2);
            writer.write(ids3);
        }

        auto lines = read_file_lines(file_path);
        REQUIRE(lines.size() == 5);
        REQUIRE(lines[0] == "FAM1\tIND1");
        REQUIRE(lines[1] == "FAM2\tIND2");
        REQUIRE(lines[2] == "FAM3\tIND3");
        REQUIRE(lines[3] == "FAM4\tIND4");
        REQUIRE(lines[4] == "FAM5\tIND5");
    }

    SECTION("Happy path - empty write then non-empty")
    {
        auto file_path = files.generate_random_file_path(".grm.id");
        std::vector<std::string> empty_ids;
        std::vector<std::string> ids = {sid("A", "B")};

        {
            GrmIdWriter writer(file_path);
            writer.write(empty_ids);
            writer.write(ids);
        }

        REQUIRE(read_file_content(file_path) == "A\tB\n");
    }
}

TEST_CASE("GrmIdWriter - Large number of IDs", "[grm_id_writer][large]")
{
    FileFixture files;
    auto file_path = files.generate_random_file_path(".grm.id");

    constexpr size_t num_ids = 1000;
    std::vector<std::string> ids;
    ids.reserve(num_ids);

    for (size_t i = 0; i < num_ids; ++i)
    {
        ids.push_back(
            sid("FAM" + std::to_string(i), "IND" + std::to_string(i)));
    }

    REQUIRE_NOTHROW(
        [&]()
        {
            GrmIdWriter writer(file_path);
            writer.write(ids);
        }());

    auto lines = read_file_lines(file_path);
    REQUIRE(lines.size() == num_ids);
    REQUIRE(lines[0] == "FAM0\tIND0");
    REQUIRE(lines[num_ids - 1] == "FAM999\tIND999");
}

TEST_CASE(
    "GrmIdWriter - IDs with special characters",
    "[grm_id_writer][special]")
{
    FileFixture files;

    SECTION("Happy path - IDs with numbers")
    {
        auto file_path = files.generate_random_file_path(".grm.id");
        std::vector<std::string> ids
            = {sid("123", "456"), sid("FAM01", "IND02")};

        {
            GrmIdWriter writer(file_path);
            writer.write(ids);
        }

        auto lines = read_file_lines(file_path);
        REQUIRE(lines.size() == 2);
        REQUIRE(lines[0] == "123\t456");
        REQUIRE(lines[1] == "FAM01\tIND02");
    }

    SECTION("Happy path - IDs with dots dashes and spaces")
    {
        auto file_path = files.generate_random_file_path(".grm.id");
        std::vector<std::string> ids
            = {sid("FAM.1", "IND-1"), sid("Project X", "Subject 42")};

        {
            GrmIdWriter writer(file_path);
            writer.write(ids);
        }

        auto lines = read_file_lines(file_path);
        REQUIRE(lines.size() == 2);
        REQUIRE(lines[0] == "FAM.1\tIND-1");
        REQUIRE(lines[1] == "Project X\tSubject 42");
    }
}
