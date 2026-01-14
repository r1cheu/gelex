#include <string_view>
#include <vector>

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_exception.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>

#include "file_fixture.h"
#include "gelex/data/sample_manager.h"
#include "gelex/exception.h"

namespace fs = std::filesystem;

using namespace gelex;  // NOLINT
using Catch::Matchers::EndsWith;
using gelex::test::FileFixture;

TEST_CASE("SampleManager - Construction and basic functionality", "[data]")
{
    FileFixture files;

    SECTION("Happy path - Construct with valid .fam file (iid_only=false)")
    {
        const auto* fam_content = R"(1 sample1 0 0 1 2.5
3 sample3 1 2 1 3.2
2 sample2 0 0 2 1.8
4 sample4 3 4 2 2.1
)";
        auto file_path = files.create_text_file(fam_content, ".fam");

        REQUIRE_NOTHROW(
            [&]()
            {
                SampleManager manager(file_path, false);

                // Verify basic properties
                REQUIRE(manager.num_common_samples() == 4);
                REQUIRE(manager.has_common_samples() == true);

                const auto& ids = manager.common_ids();
                REQUIRE(ids.size() == 4);

                // IDs should be sorted and in "FID_IID" format
                REQUIRE(ids[0] == "1_sample1");
                REQUIRE(ids[1] == "2_sample2");
                REQUIRE(ids[2] == "3_sample3");
                REQUIRE(ids[3] == "4_sample4");

                // common_id_map should be empty before finalize()
                REQUIRE(manager.common_id_map().empty());
            }());
    }

    SECTION("Happy path - Construct with valid .fam file (iid_only=true)")
    {
        const auto* fam_content = R"(1 sample1 0 0 1 2.5
3 sample3 1 2 1 3.2
2 sample2 0 0 2 1.8
4 sample4 3 4 2 2.1
)";
        auto file_path = files.create_text_file(fam_content, ".fam");

        REQUIRE_NOTHROW(
            [&]()
            {
                SampleManager manager(file_path, true);

                // Verify basic properties
                REQUIRE(manager.num_common_samples() == 4);
                REQUIRE(manager.has_common_samples() == true);

                const auto& ids = manager.common_ids();
                REQUIRE(ids.size() == 4);

                // IDs should be sorted and just IID format
                REQUIRE(ids[0] == "sample1");
                REQUIRE(ids[1] == "sample2");
                REQUIRE(ids[2] == "sample3");
                REQUIRE(ids[3] == "sample4");

                // common_id_map should be empty before finalize()
                REQUIRE(manager.common_id_map().empty());
            }());
    }
}

TEST_CASE("SampleManager - intersect() method", "[data]")
{
    FileFixture files;

    SECTION("Happy path - Intersect with overlapping IDs")
    {
        const auto* fam_content = R"(1 sample1 0 0 1 2.5
2 sample2 0 0 2 1.8
3 sample3 1 2 1 3.2
4 sample4 3 4 2 2.1
)";
        auto file_path = files.create_text_file(fam_content, ".fam");

        SampleManager manager(file_path, false);

        // Intersect with IDs that partially overlap
        std::vector<std::string> intersect_ids = {
            "2_sample2", "3_sample3", "5_sample5"  // 5_sample5 doesn't exist
        };
        manager.intersect(intersect_ids);

        // Should have only 2 common samples after intersection
        REQUIRE(manager.num_common_samples() == 2);
        REQUIRE(manager.has_common_samples() == true);

        const auto& ids = manager.common_ids();
        REQUIRE(ids.size() == 2);
        REQUIRE(ids[0] == "2_sample2");
        REQUIRE(ids[1] == "3_sample3");
    }

    SECTION("Happy path - Intersect with all matching IDs")
    {
        const auto* fam_content = R"(1 sample1 0 0 1 2.5
2 sample2 0 0 2 1.8
3 sample3 1 2 1 3.2
)";
        auto file_path = files.create_text_file(fam_content, ".fam");

        SampleManager manager(file_path, false);

        // Intersect with all existing IDs
        std::vector<std::string> intersect_ids
            = {"1_sample1", "2_sample2", "3_sample3"};
        manager.intersect(intersect_ids);

        // Should have all 3 samples
        REQUIRE(manager.num_common_samples() == 3);

        const auto& ids = manager.common_ids();
        REQUIRE(ids.size() == 3);
        REQUIRE(ids[0] == "1_sample1");
        REQUIRE(ids[1] == "2_sample2");
        REQUIRE(ids[2] == "3_sample3");
    }

    SECTION("Happy path - Intersect with no overlapping IDs")
    {
        const auto* fam_content = R"(1 sample1 0 0 1 2.5
2 sample2 0 0 2 1.8
3 sample3 1 2 1 3.2
)";
        auto file_path = files.create_text_file(fam_content, ".fam");

        SampleManager manager(file_path, false);

        // Intersect with IDs that don't exist
        std::vector<std::string> intersect_ids
            = {"4_sample4", "5_sample5"};
        manager.intersect(intersect_ids);

        // Should have no common samples after intersection
        REQUIRE(manager.num_common_samples() == 0);
        REQUIRE(manager.has_common_samples() == false);

        const auto& ids = manager.common_ids();
        REQUIRE(ids.empty());
    }

    SECTION("Edge case - Intersect with empty ID list")
    {
        const auto* fam_content = R"(1 sample1 0 0 1 2.5
2 sample2 0 0 2 1.8
)";
        auto file_path = files.create_text_file(fam_content, ".fam");

        SampleManager manager(file_path, false);

        // Initial state
        REQUIRE(manager.num_common_samples() == 2);

        // Intersect with empty list
        std::vector<std::string> intersect_ids = {};
        manager.intersect(intersect_ids);

        // Should have no samples after intersection with empty list
        REQUIRE(manager.num_common_samples() == 0);
        REQUIRE(manager.has_common_samples() == false);

        const auto& ids = manager.common_ids();
        REQUIRE(ids.empty());
    }

    SECTION("Edge case - Intersect when manager has no samples")
    {
        // Create manager with samples first
        const auto* fam_content = R"(1 sample1 0 0 1 2.5
2 sample2 0 0 2 1.8
)";
        auto file_path = files.create_text_file(fam_content, ".fam");
        SampleManager manager(file_path, false);

        // Clear samples by intersecting with empty list
        std::vector<std::string> empty_ids = {};
        manager.intersect(empty_ids);

        // Now manager has no samples
        REQUIRE(manager.num_common_samples() == 0);
        REQUIRE(manager.has_common_samples() == false);

        // Intersect with some IDs when manager has no samples
        std::vector<std::string> intersect_ids
            = {"1_sample1", "2_sample2"};
        manager.intersect(intersect_ids);

        // Should still have no samples (intersect should return early when
        // common_ids_ is empty)
        REQUIRE(manager.num_common_samples() == 0);
        REQUIRE(manager.has_common_samples() == false);
    }

    SECTION("Edge case - Intersect preserves sorting")
    {
        const auto* fam_content = R"(1 sample1 0 0 1 2.5
2 sample2 0 0 2 1.8
3 sample3 1 2 1 3.2
4 sample4 3 4 2 2.1
5 sample5 5 6 1 2.8
)";
        auto file_path = files.create_text_file(fam_content, ".fam");

        SampleManager manager(file_path, false);

        // Intersect with IDs in unsorted order
        std::vector<std::string> intersect_ids
            = {"5_sample5", "2_sample2", "4_sample4"};
        manager.intersect(intersect_ids);

        // Result should still be sorted
        REQUIRE(manager.num_common_samples() == 3);

        const auto& ids = manager.common_ids();
        REQUIRE(ids.size() == 3);
        REQUIRE(ids[0] == "2_sample2");
        REQUIRE(ids[1] == "4_sample4");
        REQUIRE(ids[2] == "5_sample5");
    }
}

TEST_CASE("SampleManager - finalize() method", "[data]")
{
    FileFixture files;

    SECTION("Happy path - finalize() creates correct mapping")
    {
        const auto* fam_content = R"(1 sample1 0 0 1 2.5
2 sample2 0 0 2 1.8
3 sample3 1 2 1 3.2
)";
        auto file_path = files.create_text_file(fam_content, ".fam");

        SampleManager manager(file_path, false);

        // Before finalize()
        REQUIRE(manager.common_id_map().empty());

        // Call finalize()
        manager.finalize();

        // After finalize()
        const auto& id_map = manager.common_id_map();
        REQUIRE(id_map.size() == 3);

        // Verify mapping is correct
        REQUIRE(id_map.at("1_sample1") == 0);
        REQUIRE(id_map.at("2_sample2") == 1);
        REQUIRE(id_map.at("3_sample3") == 2);

        // Verify consistency with common_ids()
        const auto& ids = manager.common_ids();
        for (size_t i = 0; i < ids.size(); ++i)
        {
            REQUIRE(id_map.at(ids[i]) == static_cast<Eigen::Index>(i));
        }
    }

    SECTION("Happy path - finalize() after intersect()")
    {
        const auto* fam_content = R"(1 sample1 0 0 1 2.5
2 sample2 0 0 2 1.8
3 sample3 1 2 1 3.2
4 sample4 3 4 2 2.1
)";
        auto file_path = files.create_text_file(fam_content, ".fam");

        SampleManager manager(file_path, false);

        // Intersect first
        std::vector<std::string> intersect_ids
            = {"2_sample2", "3_sample3"};
        manager.intersect(intersect_ids);

        // Then finalize
        manager.finalize();

        // Verify mapping
        const auto& id_map = manager.common_id_map();
        REQUIRE(id_map.size() == 2);
        REQUIRE(id_map.at("2_sample2") == 0);
        REQUIRE(id_map.at("3_sample3") == 1);

        // Verify consistency
        const auto& ids = manager.common_ids();
        REQUIRE(ids.size() == 2);
        REQUIRE(ids[0] == "2_sample2");
        REQUIRE(ids[1] == "3_sample3");
    }

    SECTION("Edge case - finalize() with no samples")
    {
        // Create manager with samples first
        const auto* fam_content = R"(1 sample1 0 0 1 2.5
2 sample2 0 0 2 1.8
)";
        auto file_path = files.create_text_file(fam_content, ".fam");
        SampleManager manager(file_path, false);

        // Clear samples by intersecting with empty list
        std::vector<std::string> empty_ids = {};
        manager.intersect(empty_ids);

        // Should have no samples
        REQUIRE(manager.num_common_samples() == 0);
        REQUIRE(manager.has_common_samples() == false);

        // finalize() should not crash
        REQUIRE_NOTHROW(manager.finalize());

        // Mapping should be empty
        const auto& id_map = manager.common_id_map();
        REQUIRE(id_map.empty());
    }

    SECTION("Edge case - Multiple calls to finalize()")
    {
        const auto* fam_content = R"(1 sample1 0 0 1 2.5
2 sample2 0 0 2 1.8
)";
        auto file_path = files.create_text_file(fam_content, ".fam");

        SampleManager manager(file_path, false);

        // First finalize()
        manager.finalize();

        const auto& id_map1 = manager.common_id_map();
        REQUIRE(id_map1.size() == 2);

        // Intersect and finalize again
        std::vector<std::string> intersect_ids = {"2_sample2"};
        manager.intersect(intersect_ids);
        manager.finalize();

        // Mapping should be updated
        const auto& id_map2 = manager.common_id_map();
        REQUIRE(id_map2.size() == 1);
        REQUIRE(id_map2.at("2_sample2") == 0);
    }
}

TEST_CASE("SampleManager - Integration tests", "[data]")
{
    FileFixture files;

    SECTION(
        "Happy path - Complete workflow: construct -> intersect -> finalize")
    {
        const auto* fam_content = R"(1 sample1 0 0 1 2.5
2 sample2 0 0 2 1.8
3 sample3 1 2 1 3.2
4 sample4 3 4 2 2.1
5 sample5 5 6 1 2.8
)";
        auto file_path = files.create_text_file(fam_content, ".fam");

        // 1. Construct
        SampleManager manager(file_path, false);
        REQUIRE(manager.num_common_samples() == 5);

        // 2. Intersect
        std::vector<std::string> intersect_ids = {
            "2_sample2",
            "3_sample3",
            "5_sample5",
            "6_sample6"  // 6_sample6 doesn't exist
        };
        manager.intersect(intersect_ids);
        REQUIRE(manager.num_common_samples() == 3);

        // 3. Finalize
        manager.finalize();

        // Verify final state
        REQUIRE(manager.has_common_samples() == true);

        const auto& ids = manager.common_ids();
        REQUIRE(ids.size() == 3);
        REQUIRE(ids[0] == "2_sample2");
        REQUIRE(ids[1] == "3_sample3");
        REQUIRE(ids[2] == "5_sample5");

        const auto& id_map = manager.common_id_map();
        REQUIRE(id_map.size() == 3);
        REQUIRE(id_map.at("2_sample2") == 0);
        REQUIRE(id_map.at("3_sample3") == 1);
        REQUIRE(id_map.at("5_sample5") == 2);
    }

    SECTION("Edge case - Workflow with iid_only=true")
    {
        const auto* fam_content = R"(1 sample1 0 0 1 2.5
2 sample2 0 0 2 1.8
3 sample3 1 2 1 3.2
)";
        auto file_path = files.create_text_file(fam_content, ".fam");

        // Construct with iid_only=true
        SampleManager manager(file_path, true);
        REQUIRE(manager.num_common_samples() == 3);

        // Intersect with IIDs
        std::vector<std::string> intersect_ids = {"sample2", "sample3"};
        manager.intersect(intersect_ids);
        REQUIRE(manager.num_common_samples() == 2);

        // Finalize
        manager.finalize();

        // Verify
        const auto& ids = manager.common_ids();
        REQUIRE(ids.size() == 2);
        REQUIRE(ids[0] == "sample2");
        REQUIRE(ids[1] == "sample3");

        const auto& id_map = manager.common_id_map();
        REQUIRE(id_map.size() == 2);
        REQUIRE(id_map.at("sample2") == 0);
        REQUIRE(id_map.at("sample3") == 1);
    }
}
