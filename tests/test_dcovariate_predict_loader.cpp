#include <string_view>
#include <unordered_map>
#include <vector>

#include <Eigen/Dense>
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_exception.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>

#include "../src/predictor/predict_dcovariate_loader.h"
#include "file_fixture.h"
#include "gelex/exception.h"

namespace fs = std::filesystem;

using namespace gelex::detail;  // NOLINT
using Catch::Matchers::EndsWith;
using gelex::test::FileFixture;

TEST_CASE("DcovarPredictLoader Constructor Tests", "[predictor][covariate]")
{
    FileFixture files;

    SECTION("Happy path - valid covar file with full IDs")
    {
        auto file_path = files.create_text_file(
            "FID\tIID\tSex\tPopulation\tRegion\n"
            "1\t2\tM\tEUR\tNorth\n"
            "3\t4\tF\tAFR\tSouth\n"
            "5\t6\tM\tASN\tEast\n");

        REQUIRE_NOTHROW(
            [&]()
            {
                DcovarPredictLoader loader(file_path, false);
                REQUIRE(loader.names().size() == 3);
                REQUIRE(loader.names()[0] == "Sex");
                REQUIRE(loader.names()[1] == "Population");
                REQUIRE(loader.names()[2] == "Region");

                const auto& data = loader.data();

                REQUIRE(data.size() == 3);
                REQUIRE(data.count("1_2") == 1);
                REQUIRE(data.count("3_4") == 1);
                REQUIRE(data.count("5_6") == 1);

                const auto& sample1 = data.at("1_2");
                REQUIRE(sample1.size() == 3);
                REQUIRE(sample1[0] == "M");
                REQUIRE(sample1[1] == "EUR");
                REQUIRE(sample1[2] == "North");

                const auto& sample2 = data.at("5_6");
                REQUIRE(sample2.size() == 3);
                REQUIRE(sample2[0] == "M");
                REQUIRE(sample2[1] == "ASN");
                REQUIRE(sample2[2] == "East");
            }());
    }

    SECTION("Happy path - valid covar file with IID only")
    {
        auto file_path = files.create_text_file(
            "FID\tIID\tSex\tPopulation\n"
            "1\t2\tM\tEUR\n"
            "3\t4\tF\tAFR\n");

        REQUIRE_NOTHROW(
            [&]()
            {
                DcovarPredictLoader loader(file_path, true);
                REQUIRE(loader.names().size() == 2);
                REQUIRE(loader.names()[0] == "Sex");
                REQUIRE(loader.names()[1] == "Population");
                REQUIRE(loader.data().size() == 2);

                const auto& data = loader.data();

                REQUIRE(data.size() == 2);
                REQUIRE(data.count("2") == 1);
                REQUIRE(data.count("4") == 1);

                const auto& sample1 = data.at("2");
                REQUIRE(sample1.size() == 2);
                REQUIRE(sample1[0] == "M");
                REQUIRE(sample1[1] == "EUR");

                const auto& sample2 = data.at("4");
                REQUIRE(sample2.size() == 2);
                REQUIRE(sample2[0] == "F");
                REQUIRE(sample2[1] == "AFR");
            }());
    }

    SECTION("Edge case - file with only header")
    {
        auto file_path = files.create_text_file("FID\tIID\tSex\tPopulation\n");

        REQUIRE_NOTHROW(
            [&]()
            {
                DcovarPredictLoader loader(file_path, false);
                REQUIRE(loader.names().size() == 2);
                REQUIRE(loader.data().empty());
            }());
    }

    SECTION("Exception - insufficient columns in header")
    {
        auto file_path = files.create_text_file(
            "FID\tIID\n"
            "1\t2\n");

        REQUIRE_THROWS_MATCHES(
            DcovarPredictLoader(file_path, false),
            gelex::FileFormatException,
            Catch::Matchers::MessageMatches(
                EndsWith("Covar file must have at least 3 columns, got 2")));
    }

    SECTION("Exception - column count mismatch in data row")
    {
        auto file_path = files.create_text_file(
            "FID\tIID\tSex\tPopulation\n"
            "1\t2\tM\n");  // Missing Population value

        REQUIRE_THROWS_MATCHES(
            DcovarPredictLoader(file_path, false),
            gelex::FileFormatException,
            Catch::Matchers::MessageMatches(
                EndsWith("Inconsistent number of columns at line 2")));
    }

    SECTION("Happy path - handle empty lines")
    {
        auto file_path = files.create_text_file(
            "FID\tIID\tSex\tPopulation\n"
            "\n"
            "1\t2\tM\tEUR\n"
            "\n"
            "3\t4\tF\tAFR\n"
            "\n");

        REQUIRE_NOTHROW(
            [&]()
            {
                DcovarPredictLoader loader(file_path, false);
                REQUIRE(loader.data().size() == 2);
            }());
    }

    SECTION("Edge case - single covariate")
    {
        auto file_path = files.create_text_file(
            "FID\tIID\tSex\n"
            "1\t2\tM\n"
            "3\t4\tF\n");

        REQUIRE_NOTHROW(
            [&]()
            {
                DcovarPredictLoader loader(file_path, false);
                REQUIRE(loader.names().size() == 1);
                REQUIRE(loader.names()[0] == "Sex");
                REQUIRE(loader.data().size() == 2);
            }());
    }

    SECTION("Edge case - covariate names with special characters")
    {
        auto file_path = files.create_text_file(
            "FID\tIID\tSex_Group\tPopulation-Region\n"
            "1\t2\tM\tEUR\n"
            "3\t4\tF\tAFR\n");

        REQUIRE_NOTHROW(
            [&]()
            {
                DcovarPredictLoader loader(file_path, false);
                REQUIRE(loader.names().size() == 2);
                REQUIRE(loader.names()[0] == "Sex_Group");
                REQUIRE(loader.names()[1] == "Population-Region");
                REQUIRE(loader.data().size() == 2);
            }());
    }

    SECTION("Edge case - covariate names with spaces")
    {
        auto file_path = files.create_text_file(
            "FID\tIID\tSex Group\tPopulation Region\n"
            "1\t2\tM\tEUR\n"
            "3\t4\tF\tAFR\n");

        REQUIRE_NOTHROW(
            [&]()
            {
                DcovarPredictLoader loader(file_path, false);
                REQUIRE(loader.names().size() == 2);
                REQUIRE(loader.names()[0] == "Sex Group");
                REQUIRE(loader.names()[1] == "Population Region");
                REQUIRE(loader.data().size() == 2);
            }());
    }
}

TEST_CASE("DcovarPredictLoader load Tests", "[predictor][covariate]")
{
    FileFixture files;

    SECTION("Happy path - load with complete ID mapping")
    {
        auto file_path = files.create_text_file(
            "FID\tIID\tSex\tPopulation\n"
            "1\t2\tM\tEUR\n"
            "3\t4\tF\tAFR\n"
            "5\t6\tM\tASN\n");

        DcovarPredictLoader loader(file_path, false);

        std::unordered_map<std::string, Eigen::Index> id_map
            = {{"1_2", 0}, {"3_4", 1}, {"5_6", 2}};

        auto result = loader.load(id_map);

        REQUIRE(result.size() == 2);  // Two covariates: Sex and Population
        REQUIRE(result.count("Sex") == 1);
        REQUIRE(result.count("Population") == 1);

        const auto& sex_values = result.at("Sex");
        REQUIRE(sex_values.size() == 3);
        REQUIRE(sex_values[0] == "M");
        REQUIRE(sex_values[1] == "F");
        REQUIRE(sex_values[2] == "M");

        const auto& pop_values = result.at("Population");
        REQUIRE(pop_values.size() == 3);
        REQUIRE(pop_values[0] == "EUR");
        REQUIRE(pop_values[1] == "AFR");
        REQUIRE(pop_values[2] == "ASN");
    }

    SECTION("Happy path - load with partial ID mapping")
    {
        auto file_path = files.create_text_file(
            "FID\tIID\tSex\tPopulation\n"
            "1\t2\tM\tEUR\n"
            "3\t4\tF\tAFR\n"
            "5\t6\tM\tASN\n");

        DcovarPredictLoader loader(file_path, false);

        std::unordered_map<std::string, Eigen::Index> id_map
            = {{"1_2", 0}, {"5_6", 1}};

        auto result = loader.load(id_map);

        REQUIRE(result.size() == 2);
        const auto& sex_values = result.at("Sex");
        REQUIRE(sex_values.size() == 2);
        REQUIRE(sex_values[0] == "M");  // ID "1_2" at row 0
        REQUIRE(sex_values[1] == "M");  // ID "5_6" at row 1

        const auto& pop_values = result.at("Population");
        REQUIRE(pop_values.size() == 2);
        REQUIRE(pop_values[0] == "EUR");
        REQUIRE(pop_values[1] == "ASN");
    }

    SECTION("Happy path - load with ID reordering")
    {
        auto file_path = files.create_text_file(
            "FID\tIID\tSex\tPopulation\n"
            "1\t2\tM\tEUR\n"
            "3\t4\tF\tAFR\n"
            "5\t6\tM\tASN\n");

        DcovarPredictLoader loader(file_path, false);

        std::unordered_map<std::string, Eigen::Index> id_map
            = {{"5_6", 0}, {"1_2", 1}, {"3_4", 2}};

        auto result = loader.load(id_map);

        REQUIRE(result.size() == 2);
        const auto& sex_values = result.at("Sex");
        REQUIRE(sex_values.size() == 3);
        REQUIRE(sex_values[0] == "M");  // ID "5_6" at row 0
        REQUIRE(sex_values[1] == "M");  // ID "1_2" at row 1
        REQUIRE(sex_values[2] == "F");  // ID "3_4" at row 2

        const auto& pop_values = result.at("Population");
        REQUIRE(pop_values.size() == 3);
        REQUIRE(pop_values[0] == "ASN");
        REQUIRE(pop_values[1] == "EUR");
        REQUIRE(pop_values[2] == "AFR");
    }

    SECTION("Edge case - empty ID mapping")
    {
        auto file_path = files.create_text_file(
            "FID\tIID\tSex\tPopulation\n"
            "1\t2\tM\tEUR\n");

        DcovarPredictLoader loader(file_path, false);

        std::unordered_map<std::string, Eigen::Index> id_map;

        auto result = loader.load(id_map);

        REQUIRE(result.size() == 2);  // Still two covariates
        const auto& sex_values = result.at("Sex");
        REQUIRE(sex_values.size() == 0);  // Zero samples
        const auto& pop_values = result.at("Population");
        REQUIRE(pop_values.size() == 0);
    }

    SECTION("Edge case - ID mapping with no matches")
    {
        auto file_path = files.create_text_file(
            "FID\tIID\tSex\tPopulation\n"
            "1\t2\tM\tEUR\n");

        DcovarPredictLoader loader(file_path, false);

        std::unordered_map<std::string, Eigen::Index> id_map
            = {{"nonexistent", 0}, {"another_missing", 1}};

        auto result = loader.load(id_map);

        REQUIRE(result.size() == 2);
        const auto& sex_values = result.at("Sex");
        REQUIRE(sex_values.size() == 2);
        REQUIRE(sex_values[0] == "");  // Empty string for missing ID
        REQUIRE(sex_values[1] == "");
        const auto& pop_values = result.at("Population");
        REQUIRE(pop_values.size() == 2);
        REQUIRE(pop_values[0] == "");
        REQUIRE(pop_values[1] == "");
    }

    SECTION("Happy path - load with IID only mapping")
    {
        auto file_path = files.create_text_file(
            "FID\tIID\tSex\tPopulation\n"
            "1\t2\tM\tEUR\n"
            "3\t4\tF\tAFR\n"
            "5\t6\tM\tASN\n");

        DcovarPredictLoader loader(file_path, true);

        std::unordered_map<std::string, Eigen::Index> id_map
            = {{"2", 0}, {"4", 1}, {"6", 2}};

        auto result = loader.load(id_map);

        REQUIRE(result.size() == 2);
        const auto& sex_values = result.at("Sex");
        REQUIRE(sex_values.size() == 3);
        REQUIRE(sex_values[0] == "M");
        REQUIRE(sex_values[1] == "F");
        REQUIRE(sex_values[2] == "M");

        const auto& pop_values = result.at("Population");
        REQUIRE(pop_values.size() == 3);
        REQUIRE(pop_values[0] == "EUR");
        REQUIRE(pop_values[1] == "AFR");
        REQUIRE(pop_values[2] == "ASN");
    }

    SECTION("Edge case - single covariate with multiple samples")
    {
        auto file_path = files.create_text_file(
            "FID\tIID\tGroup\n"
            "1\t2\tA\n"
            "3\t4\tB\n"
            "5\t6\tA\n");

        DcovarPredictLoader loader(file_path, false);

        std::unordered_map<std::string, Eigen::Index> id_map
            = {{"1_2", 0}, {"3_4", 1}, {"5_6", 2}};

        auto result = loader.load(id_map);

        REQUIRE(result.size() == 1);
        const auto& group_values = result.at("Group");
        REQUIRE(group_values.size() == 3);
        REQUIRE(group_values[0] == "A");
        REQUIRE(group_values[1] == "B");
        REQUIRE(group_values[2] == "A");
    }
}

TEST_CASE("DcovarPredictLoader Data Accessor Tests", "[predictor][covariate]")
{
    FileFixture files;

    SECTION("names() returns correct covariate names")
    {
        auto file_path = files.create_text_file(
            "FID\tIID\tCovar1\tCovar2\tCovar3\n"
            "1\t2\tA\tB\tC\n");

        DcovarPredictLoader loader(file_path, false);

        const auto& names = loader.names();
        REQUIRE(names.size() == 3);
        REQUIRE(names[0] == "Covar1");
        REQUIRE(names[1] == "Covar2");
        REQUIRE(names[2] == "Covar3");
    }

    SECTION("data() returns correct data mapping")
    {
        auto file_path = files.create_text_file(
            "FID\tIID\tSex\tPopulation\n"
            "1\t2\tM\tEUR\n"
            "3\t4\tF\tAFR\n");

        DcovarPredictLoader loader(file_path, false);

        const auto& data = loader.data();
        REQUIRE(data.size() == 2);
        REQUIRE(data.at("1_2") == std::vector<std::string>{"M", "EUR"});
        REQUIRE(data.at("3_4") == std::vector<std::string>{"F", "AFR"});
    }
}
