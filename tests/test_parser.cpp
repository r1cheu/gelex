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

#include <filesystem>
#include <fstream>
#include <string>
#include <vector>

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_exception.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>

#include "file_fixture.h"
#include "gelex/exception.h"
#include "gelex/io/parser.h"

namespace fs = std::filesystem;

using namespace gelex::io::detail;  // NOLINT
using Catch::Matchers::EndsWith;
using gelex::test::FileFixture;

TEST_CASE("File Stream I/O", "[parser]")
{
    FileFixture files;
    SECTION("Happy path - open existing file")
    {
        auto file_path = files.create_text_file("content");

        REQUIRE_NOTHROW(
            [&]()
            {
                auto file = open_file<std::ifstream>(file_path, std::ios::in);
                REQUIRE(file.is_open());
            }());
    }

    SECTION("Happy path - open file with custom buffer")
    {
        auto file_path = files.create_text_file("content");
        std::vector<char> buffer(4096);  // 4KB buffer

        REQUIRE_NOTHROW(
            [&]()
            {
                auto file
                    = open_file<std::ifstream>(file_path, std::ios::in, buffer);
                REQUIRE(file.is_open());
            }());
    }

    SECTION("Happy path - open file for writing")
    {
        auto file_path = files.get_test_dir() / "output.txt";

        REQUIRE_NOTHROW(
            [&]()
            {
                auto file = open_file<std::ofstream>(file_path, std::ios::out);
                REQUIRE(file.is_open());
                file << "test content";
            }());

        REQUIRE(std::filesystem::exists(file_path));
        REQUIRE(std::filesystem::file_size(file_path) > 0);
    }

    SECTION("Happy path - open file for appending")
    {
        auto file_path = files.create_text_file("initial content");

        REQUIRE_NOTHROW(
            [&]()
            {
                auto file = open_file<std::ofstream>(
                    file_path, std::ios::out | std::ios::app);
                REQUIRE(file.is_open());
                file << "\nappended content";
            }());

        // Read back to verify
        std::ifstream check(file_path);
        std::string content(
            (std::istreambuf_iterator<char>(check)),
            std::istreambuf_iterator<char>());
        REQUIRE(content.find("appended content") != std::string::npos);
    }

    SECTION("Happy path - open file for reading and writing")
    {
        auto file_path = files.create_text_file("initial content");

        REQUIRE_NOTHROW(
            [&]()
            {
                auto file = open_file<std::fstream>(
                    file_path, std::ios::in | std::ios::out);
                REQUIRE(file.is_open());

                // Read existing content
                std::string line;
                std::getline(file, line);
                REQUIRE(line == "initial content");

                // Write new content
                file.seekp(0, std::ios::end);
                file << "\nnew content";
            }());
    }

    SECTION("Exception - file not found")
    {
        REQUIRE_THROWS_MATCHES(
            open_file<std::ifstream>("non_existent_path", std::ios::in),
            gelex::GelexException,
            Catch::Matchers::MessageMatches(EndsWith("not found")));
    }

    SECTION("Exception - empty file")
    {
        auto empty_file_path = files.create_empty_file();

        REQUIRE_THROWS_MATCHES(
            open_file<std::ifstream>(empty_file_path, std::ios::in),
            gelex::GelexException,
            Catch::Matchers::MessageMatches(EndsWith("is empty")));
    }

    SECTION("Exception - directory instead of file")
    {
        auto dir_path = files.get_test_dir();

        REQUIRE_THROWS_MATCHES(
            open_file<std::ifstream>(dir_path, std::ios::in),
            gelex::GelexException,
            Catch::Matchers::MessageMatches(
                EndsWith("is a directory, not a regular file")));
        REQUIRE_THROWS_MATCHES(
            open_file<std::ofstream>(dir_path, std::ios::out),
            gelex::GelexException,
            Catch::Matchers::MessageMatches(
                EndsWith("is a directory, not a regular file")));
    }

    SECTION("Edge case - empty buffer")
    {
        auto file_path = files.create_text_file("content");
        std::vector<char> empty_buffer;

        REQUIRE_NOTHROW(
            [&]()
            {
                auto file = open_file<std::ifstream>(
                    file_path, std::ios::in, empty_buffer);
                REQUIRE(file.is_open());
            }());
    }

    SECTION("Edge case - writing to empty file is allowed")
    {
        auto empty_file_path = files.create_empty_file();

        REQUIRE_NOTHROW(
            [&]()
            {
                auto file
                    = open_file<std::ofstream>(empty_file_path, std::ios::out);
                REQUIRE(file.is_open());
                file << "writing to empty file";
            }());
    }
}

TEST_CASE("Parser Line Counting Tests", "[parser]")
{
    FileFixture files;
    SECTION("Happy path - count lines")
    {
        auto file_path = files.create_text_file("line1\nline2\nline3\nline4\n");
        REQUIRE(count_total_lines(file_path) == 4);
    }

    SECTION("Happy path - file without trailing newline")
    {
        auto file_path = files.create_text_file("line1\nline2\nline3");
        REQUIRE(count_total_lines(file_path) == 3);
    }
}
