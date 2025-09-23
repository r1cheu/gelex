#include <filesystem>
#include <fstream>
#include <iostream>
#include <vector>

#include <catch2/catch_test_macros.hpp>

#include "../src/data/loader.h"
#include "gelex/error.h"

class BimLoaderTestFixture
{
   public:
    BimLoaderTestFixture()
    {
        createValidTestFile();
        createMalformedColumnCountFile();
        createEmptyFile();
        createSingleColumnFile();
    }

    BimLoaderTestFixture(const BimLoaderTestFixture&) = default;
    BimLoaderTestFixture(BimLoaderTestFixture&&) = delete;
    BimLoaderTestFixture& operator=(const BimLoaderTestFixture&) = default;
    BimLoaderTestFixture& operator=(BimLoaderTestFixture&&) = delete;
    ~BimLoaderTestFixture()
    {
        std::remove("test_valid.bim");
        std::remove("test_malformed_columns.bim");
        std::remove("test_empty.bim");
        std::remove("test_single_column.bim");
    }

    static void createValidTestFile()
    {
        std::ofstream file("test_valid.bim");
        file << "1\trs12345\t0\t1000\tA\tG\n"
             << "1\trs67890\t0.001\t2000\tC\tT\n"
             << "2\trs24680\t0.002\t3000\tG\tA\n"
             << "X\trs13579\t0.003\t4000\tT\tC\n"
             << "1\trs11223\t0.004\t5000\tA\tT\n";
    }

    static void createMalformedColumnCountFile()
    {
        std::ofstream file("test_malformed_columns.bim");
        file << "1\trs12345\t0\t1000\tA\tG\n"
             << "1\trs67890\t0.001\t2000\tC\n";  // Missing last column
    }

    static void createEmptyFile()
    {
        std::ofstream file("test_empty.bim");
        // Empty file
    }

    static void createSingleColumnFile()
    {
        std::ofstream file("test_single_column.bim");
        file << "1\n";  // Only one column
    }
};

TEST_CASE_PERSISTENT_FIXTURE(
    BimLoaderTestFixture,
    "BimLoader::create function",
    "[loader][bim]")
{
    SECTION("Valid bim file")
    {
        auto result = gelex::detail::BimLoader::create("test_valid.bim");

        REQUIRE(result.has_value());

        const auto& snp_ids = result->ids();
        REQUIRE(snp_ids.size() == 5);
        REQUIRE(snp_ids[0] == "rs12345");
        REQUIRE(snp_ids[1] == "rs67890");
        REQUIRE(snp_ids[2] == "rs24680");
        REQUIRE(snp_ids[3] == "rs13579");
        REQUIRE(snp_ids[4] == "rs11223");
    }

    SECTION("Non-existent file")
    {
        auto result = gelex::detail::BimLoader::create("non_existent_file.bim");

        REQUIRE_FALSE(result.has_value());
        REQUIRE(result.error().code == gelex::ErrorCode::FileNotFound);
    }

    SECTION("Empty file")
    {
        auto result = gelex::detail::BimLoader::create("test_empty.bim");
        REQUIRE_FALSE(result.has_value());
        REQUIRE(result.error().code == gelex::ErrorCode::InvalidFile);
    }
}

TEST_CASE_PERSISTENT_FIXTURE(
    BimLoaderTestFixture,
    "BimLoader error handling",
    "[loader][bim]")
{
    SECTION("Malformed data - inconsistent column count")
    {
        auto result
            = gelex::detail::BimLoader::create("test_malformed_columns.bim");
        REQUIRE_FALSE(result.has_value());
        REQUIRE(result.error().code == gelex::ErrorCode::InconsistColumnCount);
    }

    SECTION("File with insufficient columns")
    {
        auto result
            = gelex::detail::BimLoader::create("test_single_column.bim");

        REQUIRE_FALSE(result.has_value());
        REQUIRE(result.error().code == gelex::ErrorCode::InconsistColumnCount);
    }
}
