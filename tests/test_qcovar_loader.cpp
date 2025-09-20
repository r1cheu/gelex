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

class QcovarLoaderTestFixture
{
   public:
    QcovarLoaderTestFixture()
    {
        // Create test files
        createValidTestFile();
        createMalformedColumnCountFile();
        createInvalidValueFile();
        createMinimalColumnsFile();
        createNoCovariatesFile();
    }

    QcovarLoaderTestFixture(const QcovarLoaderTestFixture&) = default;
    QcovarLoaderTestFixture(QcovarLoaderTestFixture&&) = delete;
    QcovarLoaderTestFixture& operator=(const QcovarLoaderTestFixture&)
        = default;
    QcovarLoaderTestFixture& operator=(QcovarLoaderTestFixture&&) = delete;
    ~QcovarLoaderTestFixture()
    {
        // Clean up test files
        std::remove("test_valid.qcovar");
        std::remove("test_malformed_columns.qcovar");
        std::remove("test_invalid_value.qcovar");
        std::remove("test_minimal.qcovar");
        std::remove("test_no_covariates.qcovar");
    }

    static void createValidTestFile()
    {
        std::ofstream file("test_valid.qcovar");
        file << "FID\tIID\tage\tbmi\theight\tweight\n"
             << "FAM1001\tIND1001\t25.5\t22.1\t175.2\t68.7\n"
             << "FAM1001\tIND1002\t30.2\t24.8\t182.5\t82.3\n"
             << "FAM1002\tIND1003\t28.7\t21.5\t168.9\t61.2\n"
             << "FAM1252\tIND1252\t35.1\t26.3\t190.1\t95.0\n";
    }

    static void createMalformedColumnCountFile()
    {
        std::ofstream file("test_malformed_columns.qcovar");
        file << "FID\tIID\tage\tbmi\theight\tweight\n"
             << "FAM1001\tIND1001\t25.5\t22.1\t175.2\t68.7\n"
             << "FAM1001\tIND1002\t30.2\t24.8\t182.5\n";  // Missing weight
                                                          // column
    }

    static void createInvalidValueFile()
    {
        std::ofstream file("test_invalid_value.qcovar");
        file << "FID\tIID\tage\tbmi\theight\tweight\n"
             << "FAM1001\tIND1001\t25.5\t22.1\t175.2\t68.7\n"
             << "FAM1001\tIND1002\t30.2\tinvalid\t182.5\t82.3\n";  // Invalid
                                                                   // bmi value
    }

    static void createMinimalColumnsFile()
    {
        std::ofstream file("test_minimal.qcovar");
        file << "FID\tIID\tage\n"
             << "FAM1001\tIND1001\t25.5\n"
             << "FAM1001\tIND1002\t30.2\n";
    }

    static void createNoCovariatesFile()
    {
        std::ofstream file("test_no_covariates.qcovar");
        file << "FID\tIID\n"
             << "FAM1001\tIND1001\n"
             << "FAM1001\tIND1002\n";
    }
};

TEST_CASE_PERSISTENT_FIXTURE(
    QcovarLoaderTestFixture,
    "QcovarLoader::create function",
    "[loader][qcovar]")
{
    SECTION("Valid qcovar file with multiple covariates")
    {
        auto result
            = gelex::detail::QcovarLoader::create("test_valid.qcovar", true);
        REQUIRE(result.has_value());

        const auto& names = result->covariate_names();
        REQUIRE(names.size() == 4);
        REQUIRE(names[0] == "age");
        REQUIRE(names[1] == "bmi");
        REQUIRE(names[2] == "height");
        REQUIRE(names[3] == "weight");

        const auto& data = result->covariate_data();
        REQUIRE(data.size() == 4);
        REQUIRE(data.contains("IND1001"));
        REQUIRE(data.contains("IND1002"));
        REQUIRE(data.contains("IND1003"));
        REQUIRE(data.contains("IND1252"));

        // Check specific values
        const auto& ind1001_data = data.at("IND1001");
        REQUIRE(ind1001_data.size() == 4);
        REQUIRE_THAT(ind1001_data[0], WithinAbs(25.5, 1e-10));
        REQUIRE_THAT(ind1001_data[1], WithinAbs(22.1, 1e-10));
        REQUIRE_THAT(ind1001_data[2], WithinAbs(175.2, 1e-10));
        REQUIRE_THAT(ind1001_data[3], WithinAbs(68.7, 1e-10));
    }

    SECTION("Valid qcovar file with minimal covariates")
    {
        auto result
            = gelex::detail::QcovarLoader::create("test_minimal.qcovar", true);

        REQUIRE(result.has_value());

        const auto& names = result->covariate_names();
        REQUIRE(names.size() == 1);
        REQUIRE(names[0] == "age");

        const auto& data = result->covariate_data();
        REQUIRE(data.size() == 2);
        REQUIRE(data.contains("IND1001"));
        REQUIRE(data.contains("IND1002"));

        const auto& ind1001_data = data.at("IND1001");
        REQUIRE(ind1001_data.size() == 1);
        REQUIRE_THAT(ind1001_data[0], WithinAbs(25.5, 1e-10));
    }

    SECTION("Non-existent file")
    {
        auto result = gelex::detail::QcovarLoader::create(
            "non_existent_file.qcovar", true);

        REQUIRE_FALSE(result.has_value());
        REQUIRE(result.error().code == gelex::ErrorCode::FileNotFound);
    }

    SECTION("IID only vs full ID mode")
    {
        auto iid_only
            = gelex::detail::QcovarLoader::create("test_valid.qcovar", true);
        REQUIRE(iid_only.has_value());
        REQUIRE(iid_only->covariate_data().contains("IND1001"));

        auto full_id
            = gelex::detail::QcovarLoader::create("test_valid.qcovar", false);
        REQUIRE(full_id.has_value());
        REQUIRE(full_id->covariate_data().contains("FAM1001_IND1001"));
    }

    SECTION("File with insufficient columns (no covariates)")
    {
        auto result = gelex::detail::QcovarLoader::create(
            "test_no_covariates.qcovar", true);

        REQUIRE_FALSE(result.has_value());
        REQUIRE(result.error().code == gelex::ErrorCode::InvalidRange);
    }
}

TEST_CASE_PERSISTENT_FIXTURE(
    QcovarLoaderTestFixture,
    "QcovarLoader::intersect method",
    "[loader][qcovar]")
{
    auto loader
        = gelex::detail::QcovarLoader::create("test_valid.qcovar", true);
    REQUIRE(loader.has_value());

    SECTION("Intersect with subset of IDs")
    {
        std::unordered_set<std::string> id_set
            = {"IND1001", "IND1002", "NON_EXISTENT", "IND1003"};

        loader->intersect(id_set);

        REQUIRE(id_set.size() == 3);
        REQUIRE(id_set.contains("IND1001"));
        REQUIRE(id_set.contains("IND1002"));
        REQUIRE(id_set.contains("IND1003"));
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
    QcovarLoaderTestFixture,
    "QcovarLoader::load method",
    "[loader][qcovar]")
{
    auto loader
        = gelex::detail::QcovarLoader::create("test_valid.qcovar", true);
    REQUIRE(loader.has_value());

    SECTION("Load with complete ID mapping")
    {
        std::unordered_map<std::string, Eigen::Index> id_map
            = {{"IND1001", 0}, {"IND1002", 2}, {"IND1003", 1}, {"IND1252", 3}};

        Eigen::MatrixXd result = loader->load(id_map);

        REQUIRE(result.rows() == 4);
        REQUIRE(result.cols() == 4);

        Eigen::MatrixXd expect{
            {25.5, 22.1, 175.2, 68.7},
            {28.7, 21.5, 168.9, 61.2},
            {30.2, 24.8, 182.5, 82.3},
            {35.1, 26.3, 190.1, 95.0}};

        REQUIRE(result.isApprox(expect, 1e-10));
    }

    SECTION("Load with partial ID mapping")
    {
        std::unordered_map<std::string, Eigen::Index> id_map
            = {{"IND1001", 1}, {"IND1252", 0}};

        Eigen::MatrixXd result = loader->load(id_map);

        REQUIRE(result.rows() == 2);
        REQUIRE(result.cols() == 4);

        Eigen::MatrixXd expect{
            {35.1, 26.3, 190.1, 95.0},  // IND1252
            {25.5, 22.1, 175.2, 68.7}   // IND1001
        };
        REQUIRE(result.isApprox(expect, 1e-10));
    }
}

TEST_CASE_PERSISTENT_FIXTURE(
    QcovarLoaderTestFixture,
    "QcovarLoader error handling",
    "[loader][qcovar]")
{
    SECTION("Malformed data - inconsistent column count")
    {
        auto result = gelex::detail::QcovarLoader::create(
            "test_malformed_columns.qcovar", true);

        REQUIRE_FALSE(result.has_value());
        REQUIRE(result.error().code == gelex::ErrorCode::InconsistColumnCount);
    }

    SECTION("Malformed data - invalid double value")
    {
        auto result = gelex::detail::QcovarLoader::create(
            "test_invalid_value.qcovar", true);

        REQUIRE_FALSE(result.has_value());
        // Should be a parsing error (NotNumber)
        REQUIRE(result.error().code == gelex::ErrorCode::NotNumber);
    }
}
