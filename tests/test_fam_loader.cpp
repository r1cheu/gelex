#include <string_view>
#include <vector>

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_exception.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>

#include "../src/data/loader/fam_loader.h"
#include "file_fixture.h"
#include "gelex/exception.h"

namespace fs = std::filesystem;

using namespace gelex::detail;  // NOLINT
using Catch::Matchers::EndsWith;
using gelex::test::FileFixture;

TEST_CASE("FamLoader - Valid .fam file loading", "[data][loader]")
{
    FileFixture files;

    SECTION(
        "Happy path - Load valid .fam file with multiple samples "
        "(iid_only=false)")
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
                FamLoader loader(file_path, false);

                const auto& ids = loader.ids();
                const auto& data = loader.data();

                REQUIRE(ids.size() == 4);
                REQUIRE(data.size() == 4);

                // Verify IDs are in "FID_IID" format
                REQUIRE(ids[0] == "1_sample1");
                REQUIRE(ids[1] == "2_sample2");
                REQUIRE(ids[2] == "3_sample3");
                REQUIRE(ids[3] == "4_sample4");

                // Verify index mapping
                REQUIRE(data.at("1_sample1") == 0);
                REQUIRE(data.at("2_sample2") == 1);
                REQUIRE(data.at("3_sample3") == 2);
                REQUIRE(data.at("4_sample4") == 3);

                for (Eigen::Index i = 0;
                     i < static_cast<Eigen::Index>(ids.size());
                     ++i)
                {
                    REQUIRE(data.at(ids[i]) == i);
                }
            }());
    }

    SECTION(
        "Happy path - Load valid .fam file with multiple samples "
        "(iid_only=true)")
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
                FamLoader loader(file_path, true);

                const auto& ids = loader.ids();
                const auto& data = loader.data();

                REQUIRE(ids.size() == 4);
                REQUIRE(data.size() == 4);

                // Verify IDs are just IID format
                REQUIRE(ids[0] == "sample1");
                REQUIRE(ids[1] == "sample2");
                REQUIRE(ids[2] == "sample3");
                REQUIRE(ids[3] == "sample4");

                // Verify index mapping
                REQUIRE(data.at("sample1") == 0);
                REQUIRE(data.at("sample2") == 1);
                REQUIRE(data.at("sample3") == 2);
                REQUIRE(data.at("sample4") == 3);
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
            FamLoader(file_path, false),
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

        FamLoader loader(file_path, false);

        // Before move
        const auto& original_ids = loader.ids();
        REQUIRE(original_ids.size() == 3);

        // Move the IDs
        auto moved_ids = std::move(loader).take_ids();
        REQUIRE(moved_ids.size() == 3);
        REQUIRE(moved_ids[0] == "1_sample1");
        REQUIRE(moved_ids[1] == "2_sample2");
        REQUIRE(moved_ids[2] == "3_sample3");

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
        REQUIRE_NOTHROW(FamLoader(file_path, true));
    }

    SECTION("Edge case - .fam file with tab delimiter (should succeed)")
    {
        const auto* fam_content = "1\tsample1\t0\t0\t1\t2.5\n";
        auto file_path = files.create_text_file(fam_content, ".fam");

        REQUIRE_NOTHROW(
            [&]()
            {
                FamLoader loader(file_path, false);
                REQUIRE(loader.ids().size() == 1);
                REQUIRE(loader.ids()[0] == "1_sample1");
            }());
    }
}
