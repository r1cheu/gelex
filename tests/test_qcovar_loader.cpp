#include <string_view>
#include <unordered_map>
#include <vector>

#include <Eigen/Dense>
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_exception.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>

#include "../src/data/loader/qcovariate_loader.h"
#include "file_fixture.h"
#include "gelex/exception.h"

namespace fs = std::filesystem;

using namespace gelex::detail;  // NOLINT
using Catch::Matchers::EndsWith;
using gelex::test::FileFixture;

TEST_CASE("QcovarLoader Constructor Tests", "[data][loader][qcovar]")
{
    FileFixture files;

    SECTION("Happy path - valid qcovar file with full IDs")
    {
        auto file_path = files.create_text_file(
            "FID\tIID\tAge\tHeight\tWeight\n"
            "1\t2\t25.5\t170.2\t65.8\n"
            "3\t4\t30.1\t165.7\t62.3\n"
            "5\t6\t28.8\t172.1\t68.9\n");

        REQUIRE_NOTHROW(
            [&]()
            {
                QuantitativeCovariateLoader loader(file_path, false);
                REQUIRE(loader.names().size() == 3);
                REQUIRE(loader.names()[0] == "Age");
                REQUIRE(loader.names()[1] == "Height");
                REQUIRE(loader.names()[2] == "Weight");

                const auto& data = loader.data();

                REQUIRE(data.size() == 3);
                REQUIRE(data.count("1_2") == 1);
                REQUIRE(data.count("3_4") == 1);
                REQUIRE(data.count("5_6") == 1);

                const auto& sample1 = data.at("1_2");
                REQUIRE(sample1.size() == 3);
                REQUIRE(sample1[0] == 25.5);
                REQUIRE(sample1[1] == 170.2);
                REQUIRE(sample1[2] == 65.8);

                const auto& sample2 = data.at("5_6");
                REQUIRE(sample2.size() == 3);
                REQUIRE(sample2[0] == 28.8);
                REQUIRE(sample2[1] == 172.1);
                REQUIRE(sample2[2] == 68.9);
            }());
    }

    SECTION("Edge case - file with only header")
    {
        auto file_path = files.create_text_file("FID\tIID\tAge\tHeight\n");

        REQUIRE_NOTHROW(
            [&]()
            {
                QuantitativeCovariateLoader loader(file_path, false);
                REQUIRE(loader.names().size() == 2);
                REQUIRE(loader.data().empty());
            }());
    }
}

TEST_CASE("QcovarLoader set_data Tests", "[data][loader][qcovar]")
{
    FileFixture files;

    SECTION("Happy path - handle empty lines")
    {
        auto file_path = files.create_text_file(
            "FID\tIID\tAge\tHeight\n"
            "\n"
            "1\t2\t25.5\t170.2\n"
            "\n"
            "3\t4\t30.1\t165.7\n"
            "\n");

        REQUIRE_NOTHROW(
            [&]()
            {
                QuantitativeCovariateLoader loader(file_path, false);
                REQUIRE(loader.data().size() == 2);
            }());
    }

    SECTION("Exception - invalid numeric data")
    {
        auto file_path = files.create_text_file(
            "FID\tIID\tAge\tHeight\n"
            "1\t2\t25.5\tinvalid\n");

        REQUIRE_THROWS_MATCHES(
            QuantitativeCovariateLoader(file_path, false),
            gelex::FileFormatException,
            Catch::Matchers::MessageMatches(
                EndsWith("failed to parse 'invalid' as number at column 3")));
    }

    SECTION("Exception - insufficient columns in data row")
    {
        auto file_path = files.create_text_file(
            "FID\tIID\tAge\tHeight\n"
            "1\t2\t25.5\n");

        REQUIRE_THROWS_MATCHES(
            QuantitativeCovariateLoader(file_path, false),
            gelex::FileFormatException,
            Catch::Matchers::MessageMatches(EndsWith(
                "expected 2 quantitative covariate values, but found 1")));
    }

    SECTION("Edge case - scientific notation numbers")
    {
        auto file_path = files.create_text_file(
            "FID\tIID\tValue1\tValue2\n"
            "1\t2\t1.23e-4\t-5.67e+3\n");

        REQUIRE_NOTHROW(
            [&]()
            {
                QuantitativeCovariateLoader loader(file_path, false);
                const auto& data = loader.data();
                REQUIRE(data.size() == 1);
                const auto& values = data.at("1_2");
                REQUIRE(values[0] == 1.23e-4);
                REQUIRE(values[1] == -5.67e+3);
            }());
    }

    SECTION("Edge case - NaN values")
    {
        auto file_path = files.create_text_file(
            "FID\tIID\tAge\tHeight\n"
            "1\t2\tnan\t170.2\n");

        REQUIRE_NOTHROW(
            [&]()
            {
                QuantitativeCovariateLoader loader(file_path, false);
                const auto& data = loader.data();
                REQUIRE(data.size() == 0);  // Row with nan should be excluded
                REQUIRE_FALSE(data.contains("1_2"));
            }());
    }

    SECTION("Edge case - Inf values")
    {
        auto file_path = files.create_text_file(
            "FID\tIID\tAge\tHeight\n"
            "1\t2\tInf\t170.2\n"
            "3\t4\t-Inf\t165.7\n"
            "5\t6\t25.5\t172.1\n");

        REQUIRE_NOTHROW(
            [&]()
            {
                QuantitativeCovariateLoader loader(file_path, false);
                const auto& data = loader.data();
                REQUIRE(data.size() == 1);
                REQUIRE(data.contains("5_6"));

                REQUIRE_FALSE(data.contains("1_2"));
                REQUIRE_FALSE(data.contains("3_4"));
                const auto& values = data.at("5_6");
                REQUIRE(values[0] == 25.5);
                REQUIRE(values[1] == 172.1);
            }());
    }
}

TEST_CASE("QcovarLoader load Tests", "[data][loader][qcovar]")
{
    FileFixture files;

    SECTION("Happy path - load with complete ID mapping")
    {
        auto file_path = files.create_text_file(
            "FID\tIID\tAge\tHeight\n"
            "1\t2\t25.5\t170.2\n"
            "3\t4\t30.1\t165.7\n"
            "5\t6\t28.8\t172.1\n");

        QuantitativeCovariateLoader loader(file_path, false);

        std::unordered_map<std::string, Eigen::Index> id_map
            = {{"1_2", 0}, {"3_4", 1}, {"5_6", 2}};

        auto qcov = loader.load(id_map);
        Eigen::MatrixXd result = std::move(qcov).X;

        REQUIRE(result.rows() == 3);
        REQUIRE(result.cols() == 2);

        REQUIRE(result(0, 0) == 25.5);
        REQUIRE(result(0, 1) == 170.2);
        REQUIRE(result(1, 0) == 30.1);
        REQUIRE(result(1, 1) == 165.7);
        REQUIRE(result(2, 0) == 28.8);
        REQUIRE(result(2, 1) == 172.1);
    }

    SECTION("Happy path - load with partial ID mapping")
    {
        auto file_path = files.create_text_file(
            "FID\tIID\tAge\tHeight\n"
            "1\t2\t25.5\t170.2\n"
            "3\t4\t30.1\t165.7\n"
            "5\t6\t28.8\t172.1\n");

        QuantitativeCovariateLoader loader(file_path, false);

        std::unordered_map<std::string, Eigen::Index> id_map = {
            {"1_2", 0}, {"5_6", 1}  // Note: "3_4" is missing from mapping
        };

        auto qcov = loader.load(id_map);
        Eigen::MatrixXd result = std::move(qcov).X;

        REQUIRE(result.rows() == 2);
        REQUIRE(result.cols() == 2);

        REQUIRE(result(0, 0) == 25.5);
        REQUIRE(result(0, 1) == 170.2);
        REQUIRE(result(1, 0) == 28.8);
        REQUIRE(result(1, 1) == 172.1);

        // Check that missing samples are not included
        REQUIRE(result.rows() == 2);
    }

    SECTION("Happy path - load with IID only mapping")
    {
        auto file_path = files.create_text_file(
            "FID\tIID\tAge\tHeight\n"
            "1\t2\t25.5\t170.2\n"
            "3\t4\t30.1\t165.7\n");

        QuantitativeCovariateLoader loader(file_path, true);

        std::unordered_map<std::string, Eigen::Index> id_map
            = {{"2", 0}, {"4", 1}};

        auto qcov = loader.load(id_map);
        Eigen::MatrixXd result = std::move(qcov).X;

        REQUIRE(result.rows() == 2);
        REQUIRE(result.cols() == 2);

        REQUIRE(result(0, 0) == 25.5);
        REQUIRE(result(0, 1) == 170.2);
        REQUIRE(result(1, 0) == 30.1);
        REQUIRE(result(1, 1) == 165.7);
    }

    SECTION("Edge case - empty ID mapping")
    {
        auto file_path = files.create_text_file(
            "FID\tIID\tAge\tHeight\n"
            "1\t2\t25.5\t170.2\n");

        QuantitativeCovariateLoader loader(file_path, false);

        std::unordered_map<std::string, Eigen::Index> id_map;

        auto qcov = loader.load(id_map);
        Eigen::MatrixXd result = std::move(qcov).X;

        REQUIRE(result.rows() == 0);
        REQUIRE(result.cols() == 2);
    }

    SECTION("Edge case - ID mapping with no matches")
    {
        auto file_path = files.create_text_file(
            "FID\tIID\tAge\tHeight\n"
            "1\t2\t25.5\t170.2\n");

        QuantitativeCovariateLoader loader(file_path, false);

        std::unordered_map<std::string, Eigen::Index> id_map
            = {{"nonexistent_id", 0}, {"another_missing", 1}};

        auto qcov = loader.load(id_map);
        Eigen::MatrixXd result = std::move(qcov).X;

        REQUIRE(result.rows() == 2);
        REQUIRE(result.cols() == 2);

        // All values should be NaN since no matches found
        REQUIRE(std::isnan(result(0, 0)));
        REQUIRE(std::isnan(result(0, 1)));
        REQUIRE(std::isnan(result(1, 0)));
        REQUIRE(std::isnan(result(1, 1)));
    }

    SECTION("Edge case - single covariate")
    {
        auto file_path = files.create_text_file(
            "FID\tIID\tAge\n"
            "1\t2\t25.5\n"
            "3\t4\t30.1\n");

        QuantitativeCovariateLoader loader(file_path, false);

        std::unordered_map<std::string, Eigen::Index> id_map
            = {{"1_2", 0}, {"3_4", 1}};

        auto qcov = loader.load(id_map);
        Eigen::MatrixXd result = std::move(qcov).X;

        REQUIRE(result.rows() == 2);
        REQUIRE(result.cols() == 1);

        REQUIRE(result(0, 0) == 25.5);
        REQUIRE(result(1, 0) == 30.1);
    }
}

TEST_CASE("QcovarLoader Integration Tests", "[data][loader][qcovar]")
{
    FileFixture files;

    SECTION("Integration - complete workflow with realistic data")
    {
        auto file_path = files.create_text_file(
            "FID\tIID\tAge\tHeight\tWeight\tBMI\n"
            "1001\t2001\t45.2\t175.3\t78.9\t25.7\n"
            "1002\t2002\t32.8\t168.4\t65.2\t23.0\n"
            "1003\t2003\t51.6\t182.1\t85.4\t25.7\n"
            "1004\t2004\t28.3\t160.9\t58.7\t22.7\n"
            "1005\t2005\t39.7\t172.8\t72.1\t24.1\n");

        QuantitativeCovariateLoader loader(file_path, false);

        // Verify names
        const auto& names = loader.names();
        REQUIRE(names.size() == 4);
        REQUIRE(names[0] == "Age");
        REQUIRE(names[1] == "Height");
        REQUIRE(names[2] == "Weight");
        REQUIRE(names[3] == "BMI");

        // Verify raw data
        const auto& data = loader.data();
        REQUIRE(data.size() == 5);

        // Create ID mapping
        std::unordered_map<std::string, Eigen::Index> id_map
            = {{"1001_2001", 0},
               {"1002_2002", 1},
               {"1003_2003", 2},
               {"1004_2004", 3},
               {"1005_2005", 4}};

        // Load into matrix
        auto qcov = loader.load(id_map);
        Eigen::MatrixXd result = std::move(qcov).X;

        REQUIRE(result.rows() == 5);
        REQUIRE(result.cols() == 4);

        // Verify specific values
        REQUIRE(result(0, 0) == 45.2);   // Sample 1 Age
        REQUIRE(result(0, 1) == 175.3);  // Sample 1 Height
        REQUIRE(result(1, 2) == 65.2);   // Sample 2 Weight
        REQUIRE(result(2, 3) == 25.7);   // Sample 3 BMI
        REQUIRE(result(3, 0) == 28.3);   // Sample 4 Age
        REQUIRE(result(4, 1) == 172.8);  // Sample 5 Height
    }

    SECTION("Integration - mixed ID formats (full and IID only)")
    {
        auto file_path = files.create_text_file(
            "FID\tIID\tCovariate1\tCovariate2\n"
            "1\t2\t1.5\t2.5\n"
            "3\t4\t3.5\t4.5\n");

        // Test with full IDs
        QuantitativeCovariateLoader loader_full(file_path, false);
        std::unordered_map<std::string, Eigen::Index> id_map_full
            = {{"1_2", 0}, {"3_4", 1}};
        auto qcov_full = loader_full.load(id_map_full);
        Eigen::MatrixXd result_full = std::move(qcov_full).X;

        // Test with IID only
        QuantitativeCovariateLoader loader_iid(file_path, true);
        std::unordered_map<std::string, Eigen::Index> id_map_iid
            = {{"2", 0}, {"4", 1}};
        auto qcov_iid = loader_iid.load(id_map_iid);
        Eigen::MatrixXd result_iid = std::move(qcov_iid).X;

        // Both should produce the same covariate values
        REQUIRE(result_full.rows() == 2);
        REQUIRE(result_iid.rows() == 2);
        REQUIRE(result_full.cols() == 2);
        REQUIRE(result_iid.cols() == 2);

        REQUIRE(result_full(0, 0) == result_iid(0, 0));
        REQUIRE(result_full(0, 1) == result_iid(0, 1));
        REQUIRE(result_full(1, 0) == result_iid(1, 0));
        REQUIRE(result_full(1, 1) == result_iid(1, 1));
    }
}
