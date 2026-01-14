#include <cmath>
#include <string_view>
#include <unordered_map>

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_exception.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>

#include "../src/data/loader/phenotype_loader.h"
#include "file_fixture.h"
#include "gelex/exception.h"

namespace fs = std::filesystem;

using namespace gelex::detail;  // NOLINT
using Catch::Matchers::EndsWith;
using gelex::test::FileFixture;

namespace
{
}  // namespace

TEST_CASE("PhenotypeLoader - Constructor Tests", "[data][phenotype]")
{
    FileFixture files;

    SECTION("Happy path - valid phenotype file with header and data")
    {
        auto file_path = files.create_text_file(
            "FID\tIID\tPhenotype1\tPhenotype2\n"
            "1\t1\t2.5\t1.0\n"
            "1\t2\t3.1\t0.5\n"
            "2\t1\t1.8\t2.2\n");

        REQUIRE_NOTHROW(
            [&]()
            {
                PhenotypeLoader loader(file_path, 2, false);
                REQUIRE(loader.name() == "Phenotype1");
                REQUIRE(loader.sample_ids().size() == 3);
                std::unordered_map<std::string, Eigen::Index> id_map
                    = {{"1_1", 0}, {"1_2", 1}, {"2_1", 2}};
                auto result = loader.load(id_map);
                REQUIRE(result(0) == 2.5);
                REQUIRE(result(1) == 3.1);
                REQUIRE(result(2) == 1.8);
            }());
        REQUIRE_NOTHROW(
            [&]()
            {
                PhenotypeLoader loader(file_path, 3, false);
                REQUIRE(loader.name() == "Phenotype2");
                REQUIRE(loader.sample_ids().size() == 3);
                std::unordered_map<std::string, Eigen::Index> id_map
                    = {{"1_1", 0}, {"1_2", 1}, {"2_1", 2}};
                auto result = loader.load(id_map);
                REQUIRE(result(0) == 1.0);
                REQUIRE(result(1) == 0.5);
                REQUIRE(result(2) == 2.2);
            }());
    }

    SECTION("Happy path - IID only mode")
    {
        auto file_path = files.create_text_file(
            "FID\tIID\tPhenotype\n"
            "1\t1\t2.5\n"
            "1\t2\t3.1\n"
            "2\t3\t1.8\n");

        REQUIRE_NOTHROW(
            [&]()
            {
                PhenotypeLoader loader(file_path, 2, true);
                REQUIRE(loader.name() == "Phenotype");
                REQUIRE(loader.sample_ids().size() == 3);
                std::unordered_map<std::string, Eigen::Index> id_map
                    = {{"1", 0}, {"2", 1}, {"3", 2}};
                auto result = loader.load(id_map);
                REQUIRE(result(0) == 2.5);
                REQUIRE(result(1) == 3.1);
                REQUIRE(result(2) == 1.8);
            }());
    }

    SECTION("Exception - column index range")
    {
        auto file_path = files.create_text_file(
            "FID\tIID\tPhenotype\n"
            "1\t1\t2.5\n");

        REQUIRE_THROWS_MATCHES(
            PhenotypeLoader(file_path, 0, false),
            gelex::FileFormatException,
            Catch::Matchers::MessageMatches(
                EndsWith("column 0 is out of range")));
        REQUIRE_THROWS_MATCHES(
            PhenotypeLoader(file_path, 5, false),
            gelex::FileFormatException,
            Catch::Matchers::MessageMatches(
                EndsWith("column 5 is out of range")));
    }
}

TEST_CASE("PhenotypeLoader - Data Loading Tests", "[data][phenotype]")
{
    FileFixture files;

    SECTION("Happy path - load with ID mapping")
    {
        auto file_path = files.create_text_file(
            "FID\tIID\tPhenotype\n"
            "1\t1\t2.5\n"
            "1\t2\t3.1\n"
            "2\t1\t1.8\n");

        PhenotypeLoader loader(file_path, 2, false);

        std::unordered_map<std::string, Eigen::Index> id_map
            = {{"1_1", 0}, {"1_2", 1}, {"2_1", 2}, {"3_1", 3}};

        auto result = loader.load(id_map);

        REQUIRE(result.size() == 4);
        REQUIRE(result(0) == 2.5);
        REQUIRE(result(1) == 3.1);
        REQUIRE(result(2) == 1.8);
        REQUIRE(std::isnan(result(3)));  // Missing phenotype
    }

    SECTION("Happy path - load with IID only mapping")
    {
        auto file_path = files.create_text_file(
            "FID\tIID\tPhenotype\n"
            "1\t1\t2.5\n"
            "1\t2\t3.1\n"
            "2\t3\t1.8\n");

        PhenotypeLoader loader(file_path, 2, true);

        std::unordered_map<std::string, Eigen::Index> id_map = {
            {"1", 0}, {"2", 1}, {"3", 2}, {"4", 3}
            // Extra ID not in phenotype data
        };

        auto result = loader.load(id_map);

        REQUIRE(result.size() == 4);
        REQUIRE(result(0) == 2.5);
        REQUIRE(result(1) == 3.1);
        REQUIRE(result(2) == 1.8);
        REQUIRE(std::isnan(result(3)));  // Missing phenotype
    }

    SECTION("Happy path - empty ID mapping")
    {
        auto file_path = files.create_text_file(
            "FID\tIID\tPhenotype\n"
            "1\t1\t2.5\n");

        PhenotypeLoader loader(file_path, 2, false);

        std::unordered_map<std::string, Eigen::Index> id_map;

        auto result = loader.load(id_map);

        REQUIRE(result.size() == 0);
    }

    SECTION("Happy path - partial ID mapping")
    {
        auto file_path = files.create_text_file(
            "FID\tIID\tPhenotype\n"
            "1\t1\t2.5\n"
            "1\t2\t3.1\n"
            "2\t1\t1.8\n");

        PhenotypeLoader loader(file_path, 2, false);

        std::unordered_map<std::string, Eigen::Index> id_map = {
            {"1_1", 0}, {"2_1", 1}  // Only two IDs from phenotype data
        };

        auto result = loader.load(id_map);

        REQUIRE(result.size() == 2);
        REQUIRE(result(0) == 2.5);
        REQUIRE(result(1) == 1.8);
    }
    SECTION("Happy path - partial ID mapping and rearrange order")
    {
        auto file_path = files.create_text_file(
            "FID\tIID\tPhenotype\n"
            "1\t1\t2.5\n"
            "1\t2\t3.1\n"
            "2\t1\t1.8\n");

        PhenotypeLoader loader(file_path, 2, false);

        std::unordered_map<std::string, Eigen::Index> id_map = {
            {"1_1", 1}, {"2_1", 0}  // Only two IDs from phenotype data
        };

        auto result = loader.load(id_map);

        REQUIRE(result.size() == 2);
        REQUIRE(result(0) == 1.8);
        REQUIRE(result(1) == 2.5);
    }
}

TEST_CASE("PhenotypeLoader - Edge Cases", "[data][phenotype]")
{
    FileFixture files;

    SECTION("Happy path - file with empty lines")
    {
        auto file_path = files.create_text_file(
            "FID\tIID\tPhenotype\n"
            "\n"  // Empty line
            "1\t1\t2.5\n"
            "\n"  // Empty line
            "1\t2\t3.1\n"
            "\n"  // Empty line
        );

        REQUIRE_NOTHROW(
            [&]()
            {
                PhenotypeLoader loader(file_path, 2, false);
                REQUIRE(loader.sample_ids().size() == 2);
                std::unordered_map<std::string, Eigen::Index> id_map
                    = {{"1_1", 0}, {"1_2", 1}};
                auto result = loader.load(id_map);
                REQUIRE(result(0) == 2.5);
                REQUIRE(result(1) == 3.1);
                REQUIRE(result.cols() == 1);
                REQUIRE(result.rows() == 2);
            }());
    }

    SECTION("Happy path - file with trailing whitespace")
    {
        auto file_path = files.create_text_file(
            "FID\tIID\tPhenotype\n"
            "1\t1\t2.5\t\n"  // Trailing tab
            "1\t2\t3.1\n"    // No trailing space
        );

        REQUIRE_NOTHROW(
            [&]()
            {
                PhenotypeLoader loader(file_path, 2, false);
                REQUIRE(loader.sample_ids().size() == 2);
                std::unordered_map<std::string, Eigen::Index> id_map
                    = {{"1_1", 0}, {"1_2", 1}};
                auto result = loader.load(id_map);
                REQUIRE(result(0) == 2.5);
                REQUIRE(result(1) == 3.1);
                REQUIRE(result.cols() == 1);
                REQUIRE(result.rows() == 2);
            }());
    }

    SECTION("Happy path - file with nan values")
    {
        auto file_path = files.create_text_file(
            "FID\tIID\tPhenotype\n"
            "1\t1\t2.5\n"
            "1\t2\tnan\n"     // nan as missing value
            "2\t1\tNaN\n"     // NaN as missing value
            "2\t2\tInf\n"     // Inf as special value
            "3\t1\t-Inf\n");  // -Inf as special value

        REQUIRE_NOTHROW(
            [&]()
            {
                PhenotypeLoader loader(file_path, 2, false);
                REQUIRE(loader.sample_ids().size() == 1);
                std::unordered_map<std::string, Eigen::Index> id_map = {
                    {"1_1", 0}, {"1_2", 1}, {"2_1", 2}, {"2_2", 3}, {"3_1", 4}};
                auto result = loader.load(id_map);
                REQUIRE(result(0) == 2.5);
                REQUIRE(std::isnan(result(1)));  // nan should be NaN in result
                REQUIRE(std::isnan(result(2)));  // NaN should be NaN in result
                REQUIRE(std::isnan(result(3)));  // Inf should be NaN in result
                REQUIRE(std::isnan(result(4)));  // -Inf should be NaN in result
            }());
    }

    SECTION("PhenotypeLoader should exclude nan values")
    {
        auto file_path = files.create_text_file(
            "FID\tIID\tPhenotype\n"
            "1\t1\t2.5\n"
            "1\t2\tnan\n"  // This should be excluded
            "2\t1\t3.1\n"
            "2\t2\tNaN\n"  // This should also be excluded
            "3\t1\t1.8\n");

        REQUIRE_NOTHROW(
            [&]()
            {
                PhenotypeLoader loader(file_path, 2, false);
                // Only non-nan values should be loaded
                REQUIRE(loader.sample_ids().size() == 3);
                std::unordered_map<std::string, Eigen::Index> id_map = {
                    {"1_1", 0}, {"1_2", 1}, {"2_1", 2}, {"2_2", 3}, {"3_1", 4}};
                auto result = loader.load(id_map);
                REQUIRE(result(0) == 2.5);
                REQUIRE(std::isnan(result(1)));  // nan should be NaN in result
                REQUIRE(result(2) == 3.1);
                REQUIRE(std::isnan(result(3)));  // NaN should be NaN in result
                REQUIRE(result(4) == 1.8);
            }());
    }

    SECTION("Exception - insufficient columns in data line")
    {
        auto file_path = files.create_text_file(
            "FID\tIID\tPhenotype\n"
            "1\t1\t2.5\n"
            "1\t2\n"  // Missing phenotype column
            "2\t1\t1.8\n");

        REQUIRE_THROWS_MATCHES(
            PhenotypeLoader(file_path, 2, false),
            gelex::FileFormatException,
            Catch::Matchers::MessageMatches(
                EndsWith("column 2 is out of range")));
    }
}
