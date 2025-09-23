#include <fstream>
#include <string>
#include <unordered_map>
#include <unordered_set>

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include "../src/data/loader.h"
#include "Eigen/Core"
#include "gelex/error.h"

using Catch::Matchers::WithinAbs;

class CovarLoaderTestFixture
{
   public:
    CovarLoaderTestFixture()
    {
        // Create test files
        createValidTestFile();
        createMalformedColumnCountFile();
        createMinimalColumnsFile();
        createNoCovariatesFile();
        createCategoricalDataFile();
    }

    CovarLoaderTestFixture(const CovarLoaderTestFixture&) = default;
    CovarLoaderTestFixture(CovarLoaderTestFixture&&) = delete;
    CovarLoaderTestFixture& operator=(const CovarLoaderTestFixture&) = default;
    CovarLoaderTestFixture& operator=(CovarLoaderTestFixture&&) = delete;
    ~CovarLoaderTestFixture()
    {
        // Clean up test files
        std::remove("test_valid.covar");
        std::remove("test_malformed_columns.covar");
        std::remove("test_minimal.covar");
        std::remove("test_no_covariates.covar");
        std::remove("test_categorical.covar");
    }

    static void createValidTestFile()
    {
        std::ofstream file("test_valid.covar");
        file << "FID\tIID\tsex\tregion\tage_group\n"
             << "FAM1001\tIND1001\tMale\tNorth\tYoung\n"
             << "FAM1001\tIND1002\tFemale\tSouth\tMiddle\n"
             << "FAM1002\tIND1003\tMale\tEast\tOld\n"
             << "FAM1252\tIND1252\tFemale\tWest\tYoung\n";
    }

    static void createCategoricalDataFile()
    {
        std::ofstream file("test_categorical.covar");
        file << "FID\tIID\tgenotype\ttreatment\tresponse\n"
             << "FAM1001\tIND1001\tAA\tControl\tGood\n"
             << "FAM1001\tIND1002\tAB\tTreatment\tExcellent\n"
             << "FAM1002\tIND1003\tBB\tControl\tFair\n"
             << "FAM1252\tIND1252\tAA\tTreatment\tGood\n"
             << "FAM1253\tIND1253\tAB\tPlacebo\tPoor\n";
    }

    static void createMalformedColumnCountFile()
    {
        std::ofstream file("test_malformed_columns.covar");
        file
            << "FID\tIID\tsex\tregion\tage_group\n"
            << "FAM1001\tIND1001\tMale\tNorth\tYoung\n"
            << "FAM1001\tIND1002\tFemale\tSouth\n";  // Missing age_group column
    }

    static void createMinimalColumnsFile()
    {
        std::ofstream file("test_minimal.covar");
        file << "FID\tIID\tsex\n"
             << "FAM1001\tIND1001\tMale\n"
             << "FAM1001\tIND1002\tFemale\n";
    }

    static void createNoCovariatesFile()
    {
        std::ofstream file("test_no_covariates.covar");
        file << "FID\tIID\n"
             << "FAM1001\tIND1001\n"
             << "FAM1001\tIND1002\n";
    }
};

TEST_CASE_PERSISTENT_FIXTURE(
    CovarLoaderTestFixture,
    "CovarLoader::create function",
    "[loader][covar]")
{
    SECTION("Valid covar file with multiple categorical covariates")
    {
        auto result
            = gelex::detail::CovarLoader::create("test_valid.covar", true);
        REQUIRE(result.has_value());

        const auto& names = result->names();
        REQUIRE(names.size() == 3);
        REQUIRE(names[0] == "sex");
        REQUIRE(names[1] == "region");
        REQUIRE(names[2] == "age_group");

        const auto& data = result->data();
        REQUIRE(data.size() == 4);
        REQUIRE(data.contains("IND1001"));
        REQUIRE(data.contains("IND1002"));
        REQUIRE(data.contains("IND1003"));
        REQUIRE(data.contains("IND1252"));

        // Check specific values
        const auto& ind1001_data = data.at("IND1001");
        REQUIRE(ind1001_data.size() == 3);
        REQUIRE(ind1001_data[0] == "Male");
        REQUIRE(ind1001_data[1] == "North");
        REQUIRE(ind1001_data[2] == "Young");
    }

    SECTION("Valid covar file with complex categorical data")
    {
        auto result = gelex::detail::CovarLoader::create(
            "test_categorical.covar", true);
        REQUIRE(result.has_value());

        const auto& names = result->names();
        REQUIRE(names.size() == 3);
        REQUIRE(names[0] == "genotype");
        REQUIRE(names[1] == "treatment");
        REQUIRE(names[2] == "response");

        const auto& data = result->data();
        REQUIRE(data.size() == 5);
    }

    SECTION("Valid covar file with minimal covariates")
    {
        auto result
            = gelex::detail::CovarLoader::create("test_minimal.covar", true);

        REQUIRE(result.has_value());

        const auto& names = result->names();
        REQUIRE(names.size() == 1);
        REQUIRE(names[0] == "sex");

        const auto& data = result->data();
        REQUIRE(data.size() == 2);
        REQUIRE(data.contains("IND1001"));
        REQUIRE(data.contains("IND1002"));

        const auto& ind1001_data = data.at("IND1001");
        REQUIRE(ind1001_data.size() == 1);
        REQUIRE(ind1001_data[0] == "Male");
    }

    SECTION("Non-existent file")
    {
        auto result = gelex::detail::CovarLoader::create(
            "non_existent_file.covar", true);

        REQUIRE_FALSE(result.has_value());
        REQUIRE(result.error().code == gelex::ErrorCode::FileNotFound);
    }

    SECTION("IID only vs full ID mode")
    {
        auto iid_only
            = gelex::detail::CovarLoader::create("test_valid.covar", true);
        REQUIRE(iid_only.has_value());
        REQUIRE(iid_only->data().contains("IND1001"));

        auto full_id
            = gelex::detail::CovarLoader::create("test_valid.covar", false);
        REQUIRE(full_id.has_value());
        REQUIRE(full_id->data().contains("FAM1001_IND1001"));
    }

    SECTION("File with insufficient columns (no covariates)")
    {
        auto result = gelex::detail::CovarLoader::create(
            "test_no_covariates.covar", true);

        REQUIRE_FALSE(result.has_value());
        REQUIRE(result.error().code == gelex::ErrorCode::InvalidRange);
    }
}

TEST_CASE_PERSISTENT_FIXTURE(
    CovarLoaderTestFixture,
    "CovarLoader::load method",
    "[loader][covar]")
{
    auto loader
        = gelex::detail::CovarLoader::create("test_categorical.covar", true);
    REQUIRE(loader.has_value());

    SECTION("Load with complete ID mapping")
    {
        std::unordered_map<std::string, Eigen::Index> id_map
            = {{"IND1001", 0},
               {"IND1002", 3},
               {"IND1003", 2},
               {"IND1252", 1},
               {"IND1253", 4}};

        Eigen::MatrixXd result = loader->load(id_map);

        REQUIRE(result.rows() == 5);
        // Should have 6 columns: genotype (2 dummy vars), treatment (2 dummy
        // vars), response (3 dummy vars)
        REQUIRE(result.cols() == 7);  // 2 + 2 + 3 = 7 dummy variables total

        Eigen::MatrixXd expected{
            {0, 0, 0, 0, 0, 1, 0},
            {0, 0, 0, 1, 0, 1, 0},
            {0, 1, 0, 0, 1, 0, 0},
            {1, 0, 0, 1, 0, 0, 0},
            {1, 0, 1, 0, 0, 0, 1},
        };

        REQUIRE(result.isApprox(expected, 1e-8));
    }

    SECTION("Load with partial ID mapping after intersection")
    {
        std::unordered_map<std::string, Eigen::Index> id_map
            = {{"IND1001", 0}, {"IND1002", 1}, {"IND1003", 2}};

        Eigen::MatrixXd result = loader->load(id_map);

        REQUIRE(result.rows() == 3);
        // After intersection,
        // genotype has 3 levels -> 2 dummy variables
        // treatment has 2 levels -> 1 dummy variable
        // response has 3 levels -> 2 dummy variables
        // Total: 2 + 1 + 2 = 5 dummy variables
        REQUIRE(result.cols() == 5);

        Eigen::MatrixXd expected{
            {0, 0, 0, 0, 1},
            {1, 0, 1, 0, 0},
            {0, 1, 0, 1, 0},
        };

        REQUIRE(result.isApprox(expected, 1e-8));
    }
}

TEST_CASE_PERSISTENT_FIXTURE(
    CovarLoaderTestFixture,
    "CovarLoader error handling",
    "[loader][covar]")
{
    SECTION("Malformed data - inconsistent column count")
    {
        auto result = gelex::detail::CovarLoader::create(
            "test_malformed_columns.covar", true);

        REQUIRE_FALSE(result.has_value());
        REQUIRE(result.error().code == gelex::ErrorCode::InconsistColumnCount);
    }
}
