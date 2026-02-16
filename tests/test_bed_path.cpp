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

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_exception.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>

#include "file_fixture.h"
#include "gelex/data/bed_path.h"
#include "gelex/exception.h"

namespace fs = std::filesystem;

using Catch::Matchers::EndsWith;
using gelex::test::FileFixture;

TEST_CASE("format_bed_path", "[data][bed_path]")
{
    FileFixture files;

    SECTION("returns existing .bed path")
    {
        auto bed_path = files.create_named_binary_file("cohort.bed", {});
        auto formatted = gelex::format_bed_path(bed_path.string());

        REQUIRE(formatted == bed_path);
    }

    SECTION("appends .bed for prefix path")
    {
        auto bed_path = files.create_named_binary_file("cohort.bed", {});
        auto prefix = bed_path;
        prefix.replace_extension();

        auto formatted = gelex::format_bed_path(prefix.string());
        REQUIRE(formatted == bed_path);
    }

    SECTION("does not confuse directory name with extension")
    {
        auto tricky_dir = files.get_test_dir() / "my.bedsets";
        fs::create_directories(tricky_dir);
        auto prefix = tricky_dir / "cohort1";
        auto bed_path = prefix;
        bed_path += ".bed";

        std::ofstream out(bed_path, std::ios::binary | std::ios::trunc);
        out.close();

        auto formatted = gelex::format_bed_path(prefix.string());
        REQUIRE(formatted == bed_path);
    }

    SECTION("throws when file does not exist")
    {
        auto missing = files.get_test_dir() / "missing_prefix";
        REQUIRE_THROWS_MATCHES(
            gelex::format_bed_path(missing.string()),
            gelex::FileNotFoundException,
            Catch::Matchers::MessageMatches(EndsWith("file not found")));
    }
}
