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

#include <string_view>
#include <vector>

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_exception.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>

#include "gelex/internal/data/loader/fam_loader.h"
#include "file_fixture.h"
#include "gelex/data/frame/dataframe_policy.h"
#include "gelex/exception.h"

namespace fs = std::filesystem;

using namespace gelex::detail;  // NOLINT
using Catch::Matchers::EndsWith;
using gelex::test::FileFixture;

namespace
{
auto sid(std::string_view fid, std::string_view iid) -> std::string
{
    return gelex::make_sample_id(fid, iid);
}
}  // namespace

TEST_CASE("FamLoader - Valid .fam file loading", "[data][loader]")
{
    FileFixture files;

    SECTION("Happy path - Load valid .fam file with multiple samples")
    {
        const auto* fam_content = R"(1 sample1 0 0 1 2.5
2 sample2 0 0 2 1.8
3 sample3 1 2 1 3.2
4 sample4 3 4 2 2.1
)";
        auto file_path = files.create_text_file(fam_content, ".fam");

        REQUIRE_NOTHROW(
            [&]()
            {
                FamLoader loader(file_path);

                const auto& ids = loader.ids();
                const auto& data = loader.data();

                REQUIRE(ids.size() == 4);
                REQUIRE(data.size() == 4);

                // Verify IDs are in canonical sample ID format
                REQUIRE(ids[0] == sid("1", "sample1"));
                REQUIRE(ids[1] == sid("2", "sample2"));
                REQUIRE(ids[2] == sid("3", "sample3"));
                REQUIRE(ids[3] == sid("4", "sample4"));

                // Verify index mapping
                REQUIRE(data.at(sid("1", "sample1")) == 0);
                REQUIRE(data.at(sid("2", "sample2")) == 1);
                REQUIRE(data.at(sid("3", "sample3")) == 2);
                REQUIRE(data.at(sid("4", "sample4")) == 3);

                for (Eigen::Index i = 0;
                     i < static_cast<Eigen::Index>(ids.size());
                     ++i)
                {
                    REQUIRE(data.at(ids[i]) == i);
                }
            }());
    }

}

TEST_CASE("FamLoader - Error conditions", "[data][loader][error]")
{
    FileFixture files;

    SECTION("Exception - Malformed .fam file (missing IID column)")
    {
        const auto* malformed_content = "1\n2 sample2\n";
        auto file_path = files.create_text_file(malformed_content, ".fam");

        REQUIRE_THROWS_MATCHES(
            FamLoader(file_path),
            gelex::FileFormatException,
            Catch::Matchers::MessageMatches(
                EndsWith("failed to parse FID and IID (missing delimiter)")));
    }
}

TEST_CASE("FamLoader - Method functionality", "[data][loader][methods]")
{
    FileFixture files;

    SECTION("Happy path - take_ids() method moves the vector")
    {
        const auto* fam_content = R"(1 sample1 0 0 1 2.5
2 sample2 0 0 2 1.8
3 sample3 1 2 1 3.2
)";
        auto file_path = files.create_text_file(fam_content, ".fam");

        FamLoader loader(file_path);

        // Before move
        const auto& original_ids = loader.ids();
        REQUIRE(original_ids.size() == 3);

        // Move the IDs
        auto moved_ids = std::move(loader).take_ids();
        REQUIRE(moved_ids.size() == 3);
        REQUIRE(moved_ids[0] == sid("1", "sample1"));
        REQUIRE(moved_ids[1] == sid("2", "sample2"));
        REQUIRE(moved_ids[2] == sid("3", "sample3"));

        // After move, the loader should still be usable but ids() may be empty
        // This depends on implementation, but at minimum shouldn't crash
        REQUIRE_NOTHROW(loader.ids());
        REQUIRE_NOTHROW(loader.data());
    }
}

TEST_CASE("FamLoader - Edge cases", "[data][loader][edge]")
{
    FileFixture files;
    SECTION("Edge case - single line")
    {
        const auto* fam_content = "1 sample1 0 0 1 2.5\n";
        auto file_path = files.create_text_file(fam_content, ".fam");
        REQUIRE_NOTHROW(FamLoader(file_path));
    }

    SECTION("Edge case - .fam file with tab delimiter (should succeed)")
    {
        const auto* fam_content = "1\tsample1\t0\t0\t1\t2.5\n";
        auto file_path = files.create_text_file(fam_content, ".fam");

        REQUIRE_NOTHROW(
            [&]()
            {
                FamLoader loader(file_path);
                REQUIRE(loader.ids().size() == 1);
                REQUIRE(loader.ids()[0] == sid("1", "sample1"));
            }());
    }
}
