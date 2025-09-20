#include <fstream>
#include <string>
#include <unordered_set>

#include <catch2/catch_test_macros.hpp>

#include "../src/data/loader.h"
#include "gelex/error.h"

class FamLoaderTestFixture
{
   public:
    FamLoaderTestFixture()
    {
        createValidTestFile();
        createMalformedColumnCountFile();
        createEmptyFile();
    }

    FamLoaderTestFixture(const FamLoaderTestFixture&) = default;
    FamLoaderTestFixture(FamLoaderTestFixture&&) = delete;
    FamLoaderTestFixture& operator=(const FamLoaderTestFixture&) = default;
    FamLoaderTestFixture& operator=(FamLoaderTestFixture&&) = delete;
    ~FamLoaderTestFixture()
    {
        std::remove("test_valid.fam");
        std::remove("test_malformed_columns.fam");
        std::remove("test_empty.fam");
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
};

TEST_CASE_PERSISTENT_FIXTURE(
    FamLoaderTestFixture,
    "FamLoader::create function",
    "[loader][fam]")
{
    SECTION("Valid fam file with IID only mode")
    {
        auto result = gelex::detail::FamLoader::create("test_valid.fam", true);

        REQUIRE(result.has_value());

        const auto& ids = result->sample_ids();
        REQUIRE(ids.size() == 5);
        REQUIRE(ids.contains("IND001"));
        REQUIRE(ids.contains("IND002"));
        REQUIRE(ids.contains("IND003"));
        REQUIRE(ids.contains("IND004"));
        REQUIRE(ids.contains("IND005"));
    }

    SECTION("Valid fam file with full ID mode")
    {
        auto result = gelex::detail::FamLoader::create("test_valid.fam", false);

        REQUIRE(result.has_value());

        const auto& ids = result->sample_ids();
        REQUIRE(ids.size() == 5);
        REQUIRE(ids.contains("FAM001_IND001"));
        REQUIRE(ids.contains("FAM001_IND002"));
        REQUIRE(ids.contains("FAM002_IND003"));
        REQUIRE(ids.contains("FAM003_IND004"));
        REQUIRE(ids.contains("FAM004_IND005"));
    }

    SECTION("Non-existent file")
    {
        auto result
            = gelex::detail::FamLoader::create("non_existent_file.fam", true);

        REQUIRE_FALSE(result.has_value());
        REQUIRE(result.error().code == gelex::ErrorCode::FileNotFound);
    }

    SECTION("Empty file")
    {
        auto result = gelex::detail::FamLoader::create("test_empty.fam", true);

        REQUIRE(result.has_value());
        REQUIRE(result->sample_ids().empty());
    }
}

TEST_CASE_PERSISTENT_FIXTURE(
    FamLoaderTestFixture,
    "FamLoader::intersect method",
    "[loader][fam]")
{
    auto loader = gelex::detail::FamLoader::create("test_valid.fam", true);
    REQUIRE(loader.has_value());

    SECTION("Intersect with subset of IDs")
    {
        std::unordered_set<std::string> id_set
            = {"IND001", "IND002", "NON_EXISTENT", "IND003"};

        loader->intersect(id_set);

        REQUIRE(id_set.size() == 3);
        REQUIRE(id_set.contains("IND001"));
        REQUIRE(id_set.contains("IND002"));
        REQUIRE(id_set.contains("IND003"));
        REQUIRE_FALSE(id_set.contains("NON_EXISTENT"));
    }

    SECTION("Intersect with empty set")
    {
        std::unordered_set<std::string> id_set;

        loader->intersect(id_set);

        REQUIRE(id_set.empty());
    }

    SECTION("Intersect with no matching IDs")
    {
        std::unordered_set<std::string> id_set
            = {"NON_EXISTENT_1", "NON_EXISTENT_2"};

        loader->intersect(id_set);

        REQUIRE(id_set.empty());
    }
}

TEST_CASE_PERSISTENT_FIXTURE(
    FamLoaderTestFixture,
    "FamLoader error handling",
    "[loader][fam]")
{
    SECTION("Malformed data - inconsistent column count")
    {
        auto result = gelex::detail::FamLoader::create(
            "test_malformed_columns.fam", true);

        REQUIRE_FALSE(result.has_value());
        REQUIRE(result.error().code == gelex::ErrorCode::InconsistColumnCount);
    }
}
