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
#include <catch2/matchers/catch_matchers_string.hpp>
#include <filesystem>
#include <fstream>
#include <iterator>
#include <string>
#include <string_view>
#include <vector>

#include "gelex/exception.h"
#include "gelex/infra/log.h"
#include "gelex/io/detail/text_writer.h"

#include "file_fixture.h"

namespace
{

namespace fs = std::filesystem;

auto read_text(const fs::path& path) -> std::string
{
    std::ifstream input{path};
    return std::string{
        std::istreambuf_iterator<char>{input},
        std::istreambuf_iterator<char>{}};
}

}  // namespace

TEST_CASE("TextWriter commits rows on destruction", "[io][text_writer]")
{
    gelex::test::FileFixture fixture;
    const auto path = fixture.get_test_dir() / "committed.tsv";
    auto temporary_path = path;
    temporary_path += ".tmp";

    {
        gelex::detail::TextWriter writer{path};
        REQUIRE(writer.path() == path);
        writer.write_header({"term", "value"});
        writer.write("alpha\t1");
        REQUIRE_FALSE(fs::exists(path));
    }

    REQUIRE(read_text(path) == "term\tvalue\nalpha\t1\n");
    REQUIRE_FALSE(fs::exists(temporary_path));
}

TEST_CASE("TextWriter discards output while unwinding", "[io][text_writer]")
{
    gelex::test::FileFixture fixture;
    const auto path = fixture.get_test_dir() / "unwound.tsv";
    auto temporary_path = path;
    temporary_path += ".tmp";

    REQUIRE_THROWS_AS(
        [&]
        {
            gelex::detail::TextWriter writer{path};
            writer.write("alpha\t1");
            throw gelex::GelexException("stop writing");
        }(),
        gelex::GelexException);

    REQUIRE_FALSE(fs::exists(path));
    REQUIRE_FALSE(fs::exists(temporary_path));
}

TEST_CASE("TextWriter reports commit failures", "[io][text_writer]")
{
    gelex::test::FileFixture fixture;
    const auto path = fixture.get_test_dir() / "failed_commit.tsv";
    auto temporary_path = path;
    temporary_path += ".tmp";
    std::vector<std::string> errors;
    gelex::set_sink(
        [&errors](gelex::Level level, std::string_view message)
        {
            if (level == gelex::Level::Error)
            {
                errors.emplace_back(message);
            }
        });
    bool obstacle_created = false;
    {
        gelex::detail::TextWriter writer{path};
        writer.write("alpha\t1");
        obstacle_created = fs::create_directory(path);
    }
    gelex::set_sink({});

    REQUIRE(obstacle_created);
    REQUIRE(fs::is_directory(path));
    REQUIRE_FALSE(fs::exists(temporary_path));
    REQUIRE(errors.size() == 1);
    CHECK_THAT(
        errors.front(), Catch::Matchers::ContainsSubstring("failed to commit"));
}
