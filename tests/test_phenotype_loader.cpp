#include <fstream>
#include <string>
#include <unordered_map>
#include <unordered_set>

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include "../src/data/loader.h"
#include "gelex/error.h"

using Catch::Matchers::WithinAbs;

class PhenotypeLoaderTestFixture
{
   public:
    PhenotypeLoaderTestFixture()
    {
        // Create test files
        createValidTestFile();
        createMalformedColumnCountFile();
        createInvalidValueFile();
    }

    PhenotypeLoaderTestFixture(const PhenotypeLoaderTestFixture&) = default;
    PhenotypeLoaderTestFixture(PhenotypeLoaderTestFixture&&) = delete;
    PhenotypeLoaderTestFixture& operator=(const PhenotypeLoaderTestFixture&)
        = default;
    PhenotypeLoaderTestFixture& operator=(PhenotypeLoaderTestFixture&&)
        = delete;
    ~PhenotypeLoaderTestFixture()
    {
        // Clean up test files
        std::remove("test_valid.phe");
        std::remove("test_malformed_columns.phe");
        std::remove("test_invalid_value.phe");
    }

    static void createValidTestFile()
    {
        std::ofstream file("test_valid.phe");
        file
            << "FID\tIID\tsex\tseason\tday\tbwt\tloc\tdam\tT1\n"
            << "FAM1001\tIND1001\tMale\tWinter\t92\t1.2\tl32\tIND0921\t4.7658\n"
            << "FAM1001\tIND1002\tMale\tSpring\t88\t2.7\tl36\tIND0921\t12."
               "4098\n"
            << "FAM1002\tIND1003\tMale\tSpring\t91\t1.0\tl17\tIND0968\t4.8545\n"
            << "FAM1252\tIND1252\tFemale\tAutumn\t82\t2.2\tl19\tIND1138\t36."
               "5418\n";
    }

    static void createMalformedColumnCountFile()
    {
        std::ofstream file("test_malformed_columns.phe");
        file
            << "FID\tIID\tsex\tseason\tday\tbwt\tloc\tdam\tT1\n"
            << "FAM1001\tIND1001\tMale\tWinter\t92\t1.2\tl32\tIND0921\t4.7658\n"
            << "FAM1001\tIND1002\tMale\tSpring\t88\t2.7\tl36\tIND0921\n";  // Missing
                                                                           // last
                                                                           // column
    }

    static void createInvalidValueFile()
    {
        std::ofstream file("test_invalid_value.phe");
        file
            << "FID\tIID\tsex\tseason\tday\tbwt\tloc\tdam\tT1\n"
            << "FAM1001\tIND1001\tMale\tWinter\t92\t1.2\tl32\tIND0921\t4.7658\n"
            << "FAM1001\tIND1002\tMale\tSpring\t88\t2.7\tl36\tIND0921\tinvalid_"
               "value\n";
    }
};

TEST_CASE_PERSISTENT_FIXTURE(
    PhenotypeLoaderTestFixture,
    "PhenotypeLoader::create function",
    "[loader][phenotype]")
{
    SECTION("Valid phenotype file with T1 column (index 8)")
    {
        auto result
            = gelex::detail::PhenotypeLoader::create("test_valid.phe", 8, true);

        REQUIRE(result.has_value());
        REQUIRE(result->name() == "T1");

        const auto& data = result->data();
        REQUIRE(data.size() == 4);

        REQUIRE_THAT(data.at("IND1001"), WithinAbs(4.7658, 1e-10));
        REQUIRE_THAT(data.at("IND1002"), WithinAbs(12.4098, 1e-10));
        REQUIRE_THAT(data.at("IND1003"), WithinAbs(4.8545, 1e-10));
        REQUIRE_THAT(data.at("IND1252"), WithinAbs(36.5418, 1e-10));
    }

    SECTION("Valid phenotype file with bwt column (index 6)")
    {
        auto result
            = gelex::detail::PhenotypeLoader::create("test_valid.phe", 5, true);

        REQUIRE(result.has_value());
        REQUIRE(result->name() == "bwt");

        const auto& data = result->data();
        REQUIRE(data.size() == 4);

        REQUIRE_THAT(data.at("IND1001"), WithinAbs(1.2, 1e-10));
        REQUIRE_THAT(data.at("IND1002"), WithinAbs(2.7, 1e-10));
        REQUIRE_THAT(data.at("IND1003"), WithinAbs(1.0, 1e-10));
        REQUIRE_THAT(data.at("IND1252"), WithinAbs(2.2, 1e-10));
    }

    SECTION("Invalid column index (too low)")
    {
        auto result
            = gelex::detail::PhenotypeLoader::create("test_valid.phe", 1, true);

        REQUIRE_FALSE(result.has_value());
        REQUIRE(result.error().code == gelex::ErrorCode::InvalidRange);
    }

    SECTION("Invalid column index (too high)")
    {
        auto result = gelex::detail::PhenotypeLoader::create(
            "test_valid.phe", 12, true);

        REQUIRE_FALSE(result.has_value());
        REQUIRE(result.error().code == gelex::ErrorCode::InvalidRange);
    }

    SECTION("Non-existent file")
    {
        auto result = gelex::detail::PhenotypeLoader::create(
            "non_existent_file.phe", 2, true);

        REQUIRE_FALSE(result.has_value());
        REQUIRE(result.error().code == gelex::ErrorCode::FileNotFound);
    }

    SECTION("IID only vs full ID mode")
    {
        auto iid_only
            = gelex::detail::PhenotypeLoader::create("test_valid.phe", 8, true);
        REQUIRE(iid_only.has_value());
        REQUIRE(iid_only->data().contains("IND1001"));

        auto full_id = gelex::detail::PhenotypeLoader::create(
            "test_valid.phe", 8, false);
        REQUIRE(full_id.has_value());
        REQUIRE(full_id->data().contains("FAM1001_IND1001"));
    }
}

TEST_CASE_PERSISTENT_FIXTURE(
    PhenotypeLoaderTestFixture,
    "PhenotypeLoader::load method",
    "[loader][phenotype]")
{
    auto loader
        = gelex::detail::PhenotypeLoader::create("test_valid.phe", 8, true);
    REQUIRE(loader.has_value());

    SECTION("Load with complete ID mapping")
    {
        std::unordered_map<std::string, Eigen::Index> id_map
            = {{"IND1001", 0}, {"IND1002", 1}, {"IND1003", 2}, {"IND1252", 3}};

        Eigen::VectorXd result = loader->load(id_map);

        REQUIRE(result.size() == 4);
        REQUIRE_THAT(result(0), WithinAbs(4.7658, 1e-10));
        REQUIRE_THAT(result(1), WithinAbs(12.4098, 1e-10));
        REQUIRE_THAT(result(2), WithinAbs(4.8545, 1e-10));
        REQUIRE_THAT(result(3), WithinAbs(36.5418, 1e-10));
    }

    SECTION("Load with disorder ID mapping")
    {
        std::unordered_map<std::string, Eigen::Index> id_map
            = {{"IND1001", 1}, {"IID1002", 0}};

        Eigen::VectorXd result = loader->load(id_map);

        REQUIRE(result.size() == 2);
        REQUIRE_THAT(result(1), WithinAbs(4.7658, 1e-10));
    }

    SECTION("Load with empty mapping")
    {
        std::unordered_map<std::string, Eigen::Index> id_map;

        Eigen::VectorXd result = loader->load(id_map);

        REQUIRE(result.size() == 0);
    }
}

TEST_CASE_PERSISTENT_FIXTURE(
    PhenotypeLoaderTestFixture,
    "PhenotypeLoader error handling",
    "[loader][phenotype]")
{
    SECTION("Malformed data - inconsistent column count")
    {
        auto result = gelex::detail::PhenotypeLoader::create(
            "test_malformed_columns.phe", 8, true);

        REQUIRE_FALSE(result.has_value());
        REQUIRE(result.error().code == gelex::ErrorCode::InconsistColumnCount);
    }

    SECTION("Malformed data - invalid double value")
    {
        auto result = gelex::detail::PhenotypeLoader::create(
            "test_invalid_value.phe", 8, true);

        REQUIRE_FALSE(result.has_value());
        // Should be a parsing error
        REQUIRE(result.error().code == gelex::ErrorCode::NotNumber);
    }
}
