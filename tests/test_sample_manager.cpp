#include <fstream>
#include <string>
#include <vector>

#include <catch2/catch_test_macros.hpp>

#include "gelex/data/sample_manager.h"
#include "gelex/error.h"

class SampleManagerTestFixture
{
   public:
    SampleManagerTestFixture()
    {
        createValidTestFile();
        createMalformedColumnCountFile();
        createEmptyFile();
        createFileWithDuplicates();
    }

    SampleManagerTestFixture(const SampleManagerTestFixture&) = default;
    SampleManagerTestFixture(SampleManagerTestFixture&&) = delete;
    SampleManagerTestFixture& operator=(const SampleManagerTestFixture&) = default;
    SampleManagerTestFixture& operator=(SampleManagerTestFixture&&) = delete;
    ~SampleManagerTestFixture()
    {
        std::remove("test_valid.fam");
        std::remove("test_malformed_columns.fam");
        std::remove("test_empty.fam");
        std::remove("test_duplicates.fam");
    }

    static void createValidTestFile()
    {
        std::ofstream file("test_valid.fam");
        file << "FAM001 IND001 0 0 1 1\n"
             << "FAM001 IND002 0 0 2 1\n"
             << "FAM002 IND003 IND001 IND002 1 2\n"
             << "FAM003 IND004 0 0 1 -9\n"
             << "FAM004 IND005 IND003 IND004 2 1\n";
    }

    static void createMalformedColumnCountFile()
    {
        std::ofstream file("test_malformed_columns.fam");
        file << "FAM001 IND001 0 0 1 1\n"
             << "FAM001 IND002 0 0 2\n";  // Missing last column
    }

    static void createEmptyFile()
    {
        std::ofstream file("test_empty.fam");
        // Empty file
    }

    static void createFileWithDuplicates()
    {
        std::ofstream file("test_duplicates.fam");
        file << "FAM001 IND001 0 0 1 1\n"
             << "FAM001 IND002 0 0 2 1\n"
             << "FAM001 IND001 0 0 1 1\n";  // Duplicate IID
    }
};

TEST_CASE_PERSISTENT_FIXTURE(
    SampleManagerTestFixture,
    "SampleManager::create function",
    "[sample_manager][create]")
{
    SECTION("Valid fam file with IID only mode")
    {
        auto result = gelex::SampleManager::create("test_valid.fam", true);

        REQUIRE(result.has_value());


        // Test accessor methods before finalize
        REQUIRE(result->num_genotyped_samples() == 5);
        REQUIRE(result->has_genotyped_samples());
        REQUIRE(result->num_common_samples() == 0);
        REQUIRE_FALSE(result->has_common_samples());

        const auto& genotyped_ids = result->genotyped_sample_ids();
        REQUIRE(genotyped_ids.size() == 5);
        REQUIRE(genotyped_ids[0] == "IND001");
        REQUIRE(genotyped_ids[1] == "IND002");
        REQUIRE(genotyped_ids[2] == "IND003");
        REQUIRE(genotyped_ids[3] == "IND004");
        REQUIRE(genotyped_ids[4] == "IND005");

        const auto& genotyped_map = result->genotyped_sample_map();
        REQUIRE(genotyped_map.size() == 5);
        REQUIRE(genotyped_map.contains("IND001"));
        REQUIRE(genotyped_map.contains("IND002"));
        REQUIRE(genotyped_map.contains("IND003"));
        REQUIRE(genotyped_map.contains("IND004"));
        REQUIRE(genotyped_map.contains("IND005"));

        // Verify indices are sequential starting from 0
        std::vector<Eigen::Index> indices;
        indices.reserve(genotyped_map.size());
        for (const auto& [id, index] : genotyped_map)
        {
            indices.push_back(index);
        }
        std::sort(indices.begin(), indices.end());
        for (size_t i = 0; i < indices.size(); ++i)
        {
            REQUIRE(indices[i] == static_cast<Eigen::Index>(i));
        }
    }

    SECTION("Valid fam file with full ID mode")
    {
        auto result = gelex::SampleManager::create("test_valid.fam", false);

        REQUIRE(result.has_value());

        // Test accessor methods before finalize
        REQUIRE(result->num_genotyped_samples() == 5);
        REQUIRE(result->has_genotyped_samples());
        REQUIRE(result->num_common_samples() == 0);
        REQUIRE_FALSE(result->has_common_samples());

        const auto& genotyped_ids = result->genotyped_sample_ids();
        REQUIRE(genotyped_ids.size() == 5);
        REQUIRE(genotyped_ids[0] == "FAM001_IND001");
        REQUIRE(genotyped_ids[1] == "FAM001_IND002");
        REQUIRE(genotyped_ids[2] == "FAM002_IND003");
        REQUIRE(genotyped_ids[3] == "FAM003_IND004");
        REQUIRE(genotyped_ids[4] == "FAM004_IND005");

        const auto& genotyped_map = result->genotyped_sample_map();
        REQUIRE(genotyped_map.size() == 5);
        REQUIRE(genotyped_map.contains("FAM001_IND001"));
        REQUIRE(genotyped_map.contains("FAM001_IND002"));
        REQUIRE(genotyped_map.contains("FAM002_IND003"));
        REQUIRE(genotyped_map.contains("FAM003_IND004"));
        REQUIRE(genotyped_map.contains("FAM004_IND005"));

        // Verify indices are sequential starting from 0
        std::vector<Eigen::Index> indices;
        indices.reserve(genotyped_map.size());
        for (const auto& [id, index] : genotyped_map)
        {
            indices.push_back(index);
        }
        std::sort(indices.begin(), indices.end());
        for (size_t i = 0; i < indices.size(); ++i)
        {
            REQUIRE(indices[i] == static_cast<Eigen::Index>(i));
        }
    }

    SECTION("Non-existent file")
    {
        auto result = gelex::SampleManager::create("non_existent_file.fam", true);

        REQUIRE_FALSE(result.has_value());
        REQUIRE(result.error().code == gelex::ErrorCode::FileNotFound);
    }

    SECTION("Empty file")
    {
        auto result = gelex::SampleManager::create("test_empty.fam", true);

        REQUIRE_FALSE(result.has_value());
        REQUIRE(result.error().code == gelex::ErrorCode::InvalidFile);
    }
}

TEST_CASE_PERSISTENT_FIXTURE(
    SampleManagerTestFixture,
    "SampleManager error handling",
    "[sample_manager][error]")
{
    SECTION("Malformed data - inconsistent column count")
    {
        auto result = gelex::SampleManager::create("test_malformed_columns.fam", true);

        REQUIRE_FALSE(result.has_value());
        REQUIRE(result.error().code == gelex::ErrorCode::InconsistColumnCount);
    }

}

TEST_CASE_PERSISTENT_FIXTURE(
    SampleManagerTestFixture,
    "SampleManager::intersect functionality",
    "[sample_manager][intersect]")
{
    SECTION("Intersect with subset of samples")
    {
        auto result = gelex::SampleManager::create("test_valid.fam", true);
        REQUIRE(result.has_value());

        std::vector<std::string_view> intersect_ids = {"IND001", "IND003", "IND005"};
        result->intersect(intersect_ids);

        // Before finalize, common_ids should be empty
        REQUIRE(result->num_common_samples() == 0);
        REQUIRE(result->num_genotyped_samples() == 5);

        result->finalize();

        // After finalize, common_ids should reflect intersection
        REQUIRE(result->num_common_samples() == 3);
        REQUIRE(result->has_common_samples());
        REQUIRE(result->num_genotyped_samples() == 5);
        REQUIRE(result->has_genotyped_samples());

        const auto& common_ids = result->common_ids();
        REQUIRE(common_ids.size() == 3);
        REQUIRE(std::ranges::find(common_ids, "IND001") != common_ids.end());
        REQUIRE(std::ranges::find(common_ids, "IND003") != common_ids.end());
        REQUIRE(std::ranges::find(common_ids, "IND005") != common_ids.end());
        REQUIRE(std::ranges::find(common_ids, "IND002") == common_ids.end());
        REQUIRE(std::ranges::find(common_ids, "IND004") == common_ids.end());

        // Verify common_ids are sorted
        REQUIRE(std::ranges::is_sorted(common_ids));

        const auto& common_map = result->common_id_map();
        REQUIRE(common_map.size() == 3);
        REQUIRE(common_map.contains("IND001"));
        REQUIRE(common_map.contains("IND003"));
        REQUIRE(common_map.contains("IND005"));

        // Verify indices are sequential starting from 0
        std::vector<Eigen::Index> indices;
        indices.reserve(common_map.size());
        for (const auto& [id, index] : common_map)
        {
            indices.push_back(index);
        }
        std::sort(indices.begin(), indices.end());
        for (size_t i = 0; i < indices.size(); ++i)
        {
            REQUIRE(indices[i] == static_cast<Eigen::Index>(i));
        }
    }

    SECTION("Intersect with empty set")
    {
        auto result = gelex::SampleManager::create("test_valid.fam", true);
        REQUIRE(result.has_value());

        std::vector<std::string_view> intersect_ids = {};
        result->intersect(intersect_ids);
        result->finalize();

        REQUIRE(result->num_common_samples() == 0);
        REQUIRE_FALSE(result->has_common_samples());
        REQUIRE(result->num_genotyped_samples() == 5);
        REQUIRE(result->has_genotyped_samples());

        const auto& common_ids = result->common_ids();
        REQUIRE(common_ids.empty());

        const auto& common_map = result->common_id_map();
        REQUIRE(common_map.empty());
    }

    SECTION("Intersect with non-matching samples")
    {
        auto result = gelex::SampleManager::create("test_valid.fam", true);
        REQUIRE(result.has_value());

        std::vector<std::string_view> intersect_ids = {"NONEXISTENT1", "NONEXISTENT2"};
        result->intersect(intersect_ids);
        result->finalize();

        REQUIRE(result->num_common_samples() == 0);
        REQUIRE_FALSE(result->has_common_samples());
        REQUIRE(result->num_genotyped_samples() == 5);
        REQUIRE(result->has_genotyped_samples());
    }

    SECTION("Intersect with partial matches")
    {
        auto result = gelex::SampleManager::create("test_valid.fam", true);
        REQUIRE(result.has_value());

        std::vector<std::string_view> intersect_ids = {"IND001", "IND003", "NONEXISTENT"};
        result->intersect(intersect_ids);
        result->finalize();

        REQUIRE(result->num_common_samples() == 2);
        REQUIRE(result->has_common_samples());

        const auto& common_ids = result->common_ids();
        REQUIRE(common_ids.size() == 2);
        REQUIRE(std::ranges::find(common_ids, "IND001") != common_ids.end());
        REQUIRE(std::ranges::find(common_ids, "IND003") != common_ids.end());
        REQUIRE(std::ranges::find(common_ids, "NONEXISTENT") == common_ids.end());
    }
}

TEST_CASE_PERSISTENT_FIXTURE(
    SampleManagerTestFixture,
    "SampleManager::finalize functionality",
    "[sample_manager][finalize]")
{
    SECTION("Finalize without intersection")
    {
        auto result = gelex::SampleManager::create("test_valid.fam", true);
        REQUIRE(result.has_value());

        // Before finalize
        REQUIRE(result->num_common_samples() == 0);
        REQUIRE(result->common_ids().empty());  // common_ids should be empty before finalize
        REQUIRE(result->common_id_map().empty());  // common_id_map should be empty before finalize

        result->finalize();

        // After finalize
        REQUIRE(result->num_common_samples() == 5);
        REQUIRE(result->has_common_samples());

        const auto& common_ids = result->common_ids();
        REQUIRE(common_ids.size() == 5);
        REQUIRE(std::ranges::find(common_ids, "IND001") != common_ids.end());
        REQUIRE(std::ranges::find(common_ids, "IND002") != common_ids.end());
        REQUIRE(std::ranges::find(common_ids, "IND003") != common_ids.end());
        REQUIRE(std::ranges::find(common_ids, "IND004") != common_ids.end());
        REQUIRE(std::ranges::find(common_ids, "IND005") != common_ids.end());

        // Verify common_ids are sorted
        REQUIRE(std::ranges::is_sorted(common_ids));

        const auto& common_map = result->common_id_map();
        REQUIRE(common_map.size() == 5);
        REQUIRE(common_map.contains("IND001"));
        REQUIRE(common_map.contains("IND002"));
        REQUIRE(common_map.contains("IND003"));
        REQUIRE(common_map.contains("IND004"));
        REQUIRE(common_map.contains("IND005"));

        // Verify indices are sequential starting from 0
        std::vector<Eigen::Index> indices;
        indices.reserve(common_map.size());
        for (const auto& [id, index] : common_map)
        {
            indices.push_back(index);
        }
        std::sort(indices.begin(), indices.end());
        for (size_t i = 0; i < indices.size(); ++i)
        {
            REQUIRE(indices[i] == static_cast<Eigen::Index>(i));
        }
    }

    SECTION("Multiple finalize calls")
    {
        auto result = gelex::SampleManager::create("test_valid.fam", true);
        REQUIRE(result.has_value());

        result->finalize();
        REQUIRE(result->num_common_samples() == 5);

        // Intersect and finalize again
        std::vector<std::string_view> intersect_ids = {"IND001", "IND003"};
        result->intersect(intersect_ids);
        result->finalize();

        REQUIRE(result->num_common_samples() == 2);
        const auto& common_ids = result->common_ids();
        REQUIRE(common_ids.size() == 2);
        REQUIRE(std::ranges::find(common_ids, "IND001") != common_ids.end());
        REQUIRE(std::ranges::find(common_ids, "IND003") != common_ids.end());
    }
}

TEST_CASE_PERSISTENT_FIXTURE(
    SampleManagerTestFixture,
    "SampleManager accessor methods",
    "[sample_manager][accessors]")
{
    SECTION("Accessors before and after operations")
    {
        auto result = gelex::SampleManager::create("test_valid.fam", true);
        REQUIRE(result.has_value());

        // Test initial state
        REQUIRE(result->num_genotyped_samples() == 5);
        REQUIRE(result->has_genotyped_samples());
        REQUIRE(result->num_common_samples() == 0);
        REQUIRE_FALSE(result->has_common_samples());

        // Test genotyped sample accessors
        const auto& genotyped_ids = result->genotyped_sample_ids();
        REQUIRE(genotyped_ids.size() == 5);
        REQUIRE(genotyped_ids[0] == "IND001");

        const auto& genotyped_map = result->genotyped_sample_map();
        REQUIRE(genotyped_map.size() == 5);
        REQUIRE(genotyped_map.contains("IND001"));

        // Test common sample accessors (should be empty before finalize)
        const auto& common_ids = result->common_ids();
        REQUIRE(common_ids.empty());

        const auto& common_map = result->common_id_map();
        REQUIRE(common_map.empty());

        // Intersect and finalize
        std::vector<std::string_view> intersect_ids = {"IND001", "IND003"};
        result->intersect(intersect_ids);
        result->finalize();

        // Test after operations
        REQUIRE(result->num_genotyped_samples() == 5);
        REQUIRE(result->has_genotyped_samples());
        REQUIRE(result->num_common_samples() == 2);
        REQUIRE(result->has_common_samples());

        // Common accessors should now be populated
        const auto& final_common_ids = result->common_ids();
        REQUIRE(final_common_ids.size() == 2);
        REQUIRE(std::ranges::find(final_common_ids, "IND001") != final_common_ids.end());
        REQUIRE(std::ranges::find(final_common_ids, "IND003") != final_common_ids.end());

        const auto& final_common_map = result->common_id_map();
        REQUIRE(final_common_map.size() == 2);
        REQUIRE(final_common_map.contains("IND001"));
        REQUIRE(final_common_map.contains("IND003"));
    }


}