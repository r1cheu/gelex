/*
 * Copyright 2026 RuLei Chen
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 * http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

#include <catch2/catch_message.hpp>
#include <catch2/matchers/catch_matchers.hpp>
#include <string_view>
#include <unordered_map>
#include <vector>

#include <Eigen/Dense>
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_exception.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>

#include "../src/data/loader/dcovariate_loader.h"
#include "file_fixture.h"
#include "gelex/exception.h"

namespace fs = std::filesystem;

using namespace gelex::detail;  // NOLINT
using Catch::Matchers::EndsWith;
using gelex::test::FileFixture;

TEST_CASE("CovarLoader Constructor Tests", "[data][loader][covar]")
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
                DiscreteCovariateLoader loader(file_path, false);
                REQUIRE(loader.sample_ids().size() == 3);
                REQUIRE(loader.sample_ids()[0] == "1_2");
                REQUIRE(loader.sample_ids()[1] == "3_4");
                REQUIRE(loader.sample_ids()[2] == "5_6");

                std::unordered_map<std::string, Eigen::Index> id_map
                    = {{"1_2", 0}, {"3_4", 1}, {"5_6", 2}};

                auto dcov = loader.load(id_map);
                // Sex: F (baseline), M -> 1
                // Population: AFR (baseline), ASN, EUR -> 2
                // Region: East (baseline), North, South -> 2
                // Total columns: 5
                REQUIRE(dcov.X.cols() == 5);
                REQUIRE(dcov.X.rows() == 3);

                // Sample 1: 1_2 (M, EUR, North)
                // Sex_M=1, Pop_ASN=0, Pop_EUR=1, Reg_North=1, Reg_South=0
                REQUIRE(dcov.X(0, 0) == 1.0);
                REQUIRE(dcov.X(0, 1) == 0.0);
                REQUIRE(dcov.X(0, 2) == 1.0);
                REQUIRE(dcov.X(0, 3) == 1.0);
                REQUIRE(dcov.X(0, 4) == 0.0);
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
                DiscreteCovariateLoader loader(file_path, true);
                REQUIRE(loader.sample_ids().size() == 2);
                REQUIRE(loader.sample_ids()[0] == "2");
                REQUIRE(loader.sample_ids()[1] == "4");

                std::unordered_map<std::string, Eigen::Index> id_map
                    = {{"2", 0}, {"4", 1}};
                auto dcov = loader.load(id_map);

                // Sex: F (baseline), M -> 1
                // Population: AFR (baseline), EUR -> 1
                REQUIRE(dcov.X.cols() == 2);
                REQUIRE(dcov.X.rows() == 2);

                // Sample 1: 2 (M, EUR) -> [1, 1]
                REQUIRE(dcov.X(0, 0) == 1.0);
                REQUIRE(dcov.X(0, 1) == 1.0);
            }());
    }

    SECTION("Edge case - file with only header")
    {
        auto file_path = files.create_text_file("FID\tIID\tSex\tPopulation\n");

        REQUIRE_NOTHROW(
            [&]()
            {
                DiscreteCovariateLoader loader(file_path, false);
                REQUIRE(loader.sample_ids().empty());
            }());
    }

    SECTION("Exception - insufficient columns in header")
    {
        auto file_path = files.create_text_file(
            "FID\tIID\n"
            "1\t2\n");

        REQUIRE_THROWS_MATCHES(
            DiscreteCovariateLoader(file_path, false),
            gelex::FileFormatException,
            Catch::Matchers::MessageMatches(
                EndsWith("categorical covariates must have > 2 columns")));
    }
}

TEST_CASE("CovarLoader set_data Tests", "[data][loader][covar]")
{
    FileFixture files;

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
                DiscreteCovariateLoader loader(file_path, false);
                REQUIRE(loader.sample_ids().size() == 2);
            }());
    }

    SECTION("Happy path - handle missing categorical values")
    {
        auto file_path = files.create_text_file(
            "FID\tIID\tSex\tPopulation\tRegion\n"
            "1\t2\tM\tEUR\t\n"
            "3\t4\tF\t\tSouth\n"
            "5\t6\t\tASN\tEast\n");

        REQUIRE_THROWS_MATCHES(
            DiscreteCovariateLoader(file_path, false),
            gelex::FileFormatException,
            Catch::Matchers::MessageMatches(
                EndsWith("empty value encountered")));
    }

    SECTION("Exception - column count mismatch in data row")
    {
        auto file_path = files.create_text_file(
            "FID\tIID\tSex\tPopulation\n"
            "1\t2\tM\n");  // Missing Population value

        REQUIRE_THROWS_MATCHES(
            DiscreteCovariateLoader(file_path, false),
            gelex::FileFormatException,
            Catch::Matchers::MessageMatches(EndsWith("Column count mismatch")));
    }

    SECTION("Edge case - single categorical covariate")
    {
        auto file_path = files.create_text_file(
            "FID\tIID\tSex\n"
            "1\t2\tM\n"
            "3\t4\tF\n");

        REQUIRE_NOTHROW(
            [&]()
            {
                DiscreteCovariateLoader loader(file_path, false);
                REQUIRE(loader.sample_ids().size() == 2);
                auto dcov = loader.load({{"1_2", 0}, {"3_4", 1}});
                REQUIRE(dcov.X.cols() == 1);
                REQUIRE(dcov.X(0, 0) == 1);
            }());
    }
}

TEST_CASE("CovarLoader load Tests", "[data][loader][covar]")
{
    FileFixture files;

    SECTION("Happy path - load with complete ID mapping and one-hot encoding")
    {
        auto file_path = files.create_text_file(
            "FID\tIID\tSex\tPopulation\n"
            "1\t2\tM\tEUR\n"
            "3\t4\tF\tAFR\n"
            "5\t6\tM\tASN\n");

        DiscreteCovariateLoader loader(file_path, false);

        std::unordered_map<std::string, Eigen::Index> id_map
            = {{"1_2", 0}, {"3_4", 1}, {"5_6", 2}};

        auto dcov = loader.load(id_map);
        Eigen::MatrixXd result = std::move(dcov).X;

        // Expected encoding:
        // Sex: M (baseline), F -> 1 dummy variable
        // Population: AFR(baseline), ASN, EUR-> 2 dummy variables
        // Total columns: 1 + 2 = 3
        REQUIRE(result.rows() == 3);
        REQUIRE(result.cols() == 3);

        // Sample 1: M, EUR -> [1, 0, 1] (all baseline)
        REQUIRE(result(0, 0) == 1.0);  // Sex_F
        REQUIRE(result(0, 1) == 0.0);  // Population_ASN
        REQUIRE(result(0, 2) == 1.0);  // Population_EUR

        // Sample 2: F, AFR -> [0, 0, 0]
        REQUIRE(result(1, 0) == 0.0);  // Sex_F
        REQUIRE(result(1, 1) == 0.0);  // Population_ASN
        REQUIRE(result(1, 2) == 0.0);  // Population_EUR

        // Sample 3: M, ASN -> [1, 1, 0]
        REQUIRE(result(2, 0) == 1.0);  // Sex_F
        REQUIRE(result(2, 1) == 1.0);  // Population_ASN
        REQUIRE(result(2, 2) == 0.0);  // Population_EUR
    }

    SECTION("Happy path - load with partial ID mapping")
    {
        auto file_path = files.create_text_file(
            "FID\tIID\tSex\tPopulation\n"
            "1\t2\tM\tEUR\n"
            "3\t4\tF\tAFR\n"
            "5\t6\tM\tASN\n");

        DiscreteCovariateLoader loader(file_path, false);

        std::unordered_map<std::string, Eigen::Index> id_map
            = {{"1_2", 0}, {"5_6", 1}};

        auto dcov = loader.load(id_map);
        Eigen::MatrixXd result = std::move(dcov).X;

        REQUIRE(result.rows() == 2);
        REQUIRE(result.cols() == 1);

        // Sample 1: M, EUR -> [1]
        REQUIRE(result(0, 0) == 1);  // Population_ASN
        // Sample 3: M, ASN -> [0]
        REQUIRE(result(1, 0) == 0);  // Population_ASN
    }
    SECTION("Happy path - load with partial ID mapping and reordering")
    {
        auto file_path = files.create_text_file(
            "FID\tIID\tSex\tPopulation\n"
            "1\t2\tM\tEUR\n"
            "3\t4\tF\tAFR\n"
            "5\t6\tM\tASN\n");

        DiscreteCovariateLoader loader(file_path, false);

        std::unordered_map<std::string, Eigen::Index> id_map
            = {{"1_2", 1}, {"5_6", 0}};

        auto dcov = loader.load(id_map);
        Eigen::MatrixXd result = std::move(dcov).X;

        REQUIRE(result.rows() == 2);
        REQUIRE(result.cols() == 1);

        // Sample 3: M, ASN -> [0]
        REQUIRE(result(0, 0) == 0);  // Population_ASN
        // Sample 1: M, EUR -> [1]
        REQUIRE(result(1, 0) == 1);  // Population_ASN
    }

    SECTION("Happy path - load with IID only mapping")
    {
        auto file_path = files.create_text_file(
            "FID\tIID\tSex\tPopulation\n"
            "1\t2\tM\tEUR\n"
            "3\t4\tF\tAFR\n"
            "5\t6\tM\tASN\n");

        DiscreteCovariateLoader loader(file_path, true);

        std::unordered_map<std::string, Eigen::Index> id_map
            = {{"2", 0}, {"4", 1}, {"6", 2}};

        auto dcov = loader.load(id_map);
        Eigen::MatrixXd result = std::move(dcov).X;

        // Expected encoding:
        // Sex: M (baseline), F -> 1 dummy variable
        // Population: AFR(baseline), ASN, EUR-> 2 dummy variables
        // Total columns: 1 + 2 = 3
        REQUIRE(result.rows() == 3);
        REQUIRE(result.cols() == 3);

        // Sample 1: M, EUR -> [1, 0, 1] (all baseline)
        REQUIRE(result(0, 0) == 1.0);  // Sex_F
        REQUIRE(result(0, 1) == 0.0);  // Population_ASN
        REQUIRE(result(0, 2) == 1.0);  // Population_EUR

        // Sample 2: F, AFR -> [0, 0, 0]
        REQUIRE(result(1, 0) == 0.0);  // Sex_F
        REQUIRE(result(1, 1) == 0.0);  // Population_ASN
        REQUIRE(result(1, 2) == 0.0);  // Population_EUR

        // Sample 3: M, ASN -> [1, 1, 0]
        REQUIRE(result(2, 0) == 1.0);  // Sex_F
        REQUIRE(result(2, 1) == 1.0);  // Population_ASN
        REQUIRE(result(2, 2) == 0.0);  // Population_EUR
    }

    SECTION("Edge case - empty ID mapping")
    {
        auto file_path = files.create_text_file(
            "FID\tIID\tSex\tPopulation\n"
            "1\t2\tM\tEUR\n");

        DiscreteCovariateLoader loader(file_path, false);

        std::unordered_map<std::string, Eigen::Index> id_map;

        auto dcov = loader.load(id_map);
        Eigen::MatrixXd result = std::move(dcov).X;

        REQUIRE(result.rows() == 0);
        REQUIRE(result.cols() == 0);  // No valid samples, so no dummy variables
    }

    SECTION("Edge case - ID mapping with no matches")
    {
        auto file_path = files.create_text_file(
            "FID\tIID\tSex\tPopulation\n"
            "1\t2\tM\tEUR\n");

        DiscreteCovariateLoader loader(file_path, false);

        std::unordered_map<std::string, Eigen::Index> id_map
            = {{"nonexistent_id", 0}, {"another_missing", 1}};

        auto dcov = loader.load(id_map);
        Eigen::MatrixXd result = std::move(dcov).X;

        REQUIRE(result.rows() == 2);
        REQUIRE(result.cols() == 0);
    }

    SECTION("Edge case - single categorical variable with two levels")
    {
        auto file_path = files.create_text_file(
            "FID\tIID\tSex\n"
            "1\t2\tM\n"
            "3\t4\tF\n");

        DiscreteCovariateLoader loader(file_path, false);

        std::unordered_map<std::string, Eigen::Index> id_map
            = {{"1_2", 0}, {"3_4", 1}};

        auto dcov = loader.load(id_map);
        Eigen::MatrixXd result = std::move(dcov).X;

        REQUIRE(result.rows() == 2);
        REQUIRE(result.cols() == 1);  // One dummy variable for Sex

        REQUIRE(result(0, 0) == 1);  // M
        REQUIRE(result(1, 0) == 0);  // F (baseline)
    }

    SECTION("Edge case - categorical variable with single level")
    {
        auto file_path = files.create_text_file(
            "FID\tIID\tSex\n"
            "1\t2\tM\n"
            "3\t4\tM\n");  // Only one level

        DiscreteCovariateLoader loader(file_path, false);

        std::unordered_map<std::string, Eigen::Index> id_map
            = {{"1_2", 0}, {"3_4", 1}};

        auto dcov = loader.load(id_map);
        Eigen::MatrixXd result = std::move(dcov).X;

        REQUIRE(result.rows() == 2);
        REQUIRE(result.cols() == 0);  // No dummy variables for single level
        REQUIRE(result.isZero());     // All zeros
    }

    SECTION("Edge case - categorical variable with missing values")
    {
        auto file_path = files.create_text_file(
            "FID\tIID\tSex\tPopulation\n"
            "1\t2\tM\tEUR\n"
            "3\t4\t\tAFR\n"  // Missing Sex
            "5\t6\tF\t\n");  // Missing Population

        REQUIRE_THROWS_MATCHES(
            DiscreteCovariateLoader(file_path, false),
            gelex::FileFormatException,
            Catch::Matchers::MessageMatches(
                EndsWith("empty value encountered")));
    }
}

TEST_CASE("CovarLoader Integration Tests", "[data][loader][covar]")
{
    FileFixture files;
    SECTION(
        "Integration - complex categorical encoding with numeric-like values")
    {
        auto file_path = files.create_text_file(
            "FID\tIID\tGroup\tCategory\n"
            "1\t2\t1\tA\n"
            "3\t4\t2\tB\n"
            "5\t6\t1\tC\n"
            "7\t8\t3\tA\n");

        DiscreteCovariateLoader loader(file_path, false);

        std::unordered_map<std::string, Eigen::Index> id_map
            = {{"1_2", 0}, {"3_4", 1}, {"5_6", 2}, {"7_8", 3}};

        auto dcov = loader.load(id_map);
        Eigen::MatrixXd result = std::move(dcov).X;

        REQUIRE(result.rows() == 4);
        // Group: 1 (baseline), 2, 3 -> 2 dummies
        // Category: A (baseline), B, C -> 2 dummies
        // Total: 2 + 2 = 4 columns
        REQUIRE(result.cols() == 4);

        // Sample 1: Group=1, Category=A -> [0,0,0,0]
        REQUIRE(result.row(0).isZero());

        // Sample 2: Group=2, Category=B -> [1,0,1,0]
        REQUIRE(result(1, 0) == 1.0);  // Group_2
        REQUIRE(result(1, 1) == 0.0);  // Group_3
        REQUIRE(result(1, 2) == 1.0);  // Category_B
        REQUIRE(result(1, 3) == 0.0);  // Category_C

        // Sample 3: Group=1, Category=C -> [0,0,0,1]
        REQUIRE(result(2, 0) == 0.0);  // Group_2
        REQUIRE(result(2, 1) == 0.0);  // Group_3
        REQUIRE(result(2, 2) == 0.0);  // Category_B
        REQUIRE(result(2, 3) == 1.0);  // Category_C

        // Sample 4: Group=3, Category=A -> [0,1,0,0]
        REQUIRE(result(3, 0) == 0.0);  // Group_2
        REQUIRE(result(3, 1) == 1.0);  // Group_3
        REQUIRE(result(3, 2) == 0.0);  // Category_B
        REQUIRE(result(3, 3) == 0.0);  // Category_C
    }
}

TEST_CASE("CovarLoader nan/inf exclusion tests", "[data][loader][covar]")
{
    FileFixture files;

    SECTION("Edge case - exclude rows with nan/inf string values")
    {
        auto file_path = files.create_text_file(
            "FID\tIID\tSex\tGroup\tCategory\n"
            "1\t2\tM\tnan\tA\n"
            "3\t4\tF\tB\tNaN\n"
            "5\t6\tM\tC\tinf\n"
            "7\t8\tF\tD\tInf\n"
            "9\t10\tM\tE\t-inf\n"
            "11\t12\tF\tF\t+Inf\n"
            "13\t14\tM\tG\tValid\n");

        REQUIRE_NOTHROW(
            [&]()
            {
                DiscreteCovariateLoader loader(file_path, false);
                const auto& ids = loader.sample_ids();
                // only the last row should be retained
                REQUIRE(ids.size() == 1);
                REQUIRE(ids[0] == "13_14");

                auto dcov = loader.load({{"13_14", 0}});
                REQUIRE(dcov.X.cols() == 0);  // All single levels
            }());
    }

    SECTION("Edge case - mixed valid and invalid values in row")
    {
        auto file_path = files.create_text_file(
            "FID\tIID\tSex\tGroup\n"
            "1\t2\tM\tnan\n"
            "3\t4\tF\tB\n"
            "5\t6\tinf\tC\n");

        REQUIRE_NOTHROW(
            [&]()
            {
                DiscreteCovariateLoader loader(file_path, false);
                const auto& ids = loader.sample_ids();
                REQUIRE(ids.size() == 1);
                REQUIRE(ids[0] == "3_4");
            }());
    }

    SECTION("Integration - nan/inf values not included in encoding levels")
    {
        auto file_path = files.create_text_file(
            "FID\tIID\tGroup\n"
            "1\t2\tnan\n"
            "3\t4\tA\n"
            "5\t6\tinf\n"
            "7\t8\tA\n");

        DiscreteCovariateLoader loader(file_path, false);
        std::unordered_map<std::string, Eigen::Index> id_map
            = {{"3_4", 0}, {"7_8", 1}};
        auto dcov = loader.load(id_map);
        Eigen::MatrixXd result = std::move(dcov).X;

        // 只有 "A" 一个水平，没有 dummy 变量
        REQUIRE(result.rows() == 2);
        REQUIRE(result.cols() == 0);
    }
}
