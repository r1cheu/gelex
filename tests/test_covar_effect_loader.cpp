#include <string_view>
#include <unordered_map>

#include <Eigen/Dense>
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_exception.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>

#include "../src/predictor/covar_effect_loader.h"
#include "file_fixture.h"
#include "gelex/exception.h"

namespace fs = std::filesystem;

using namespace gelex::detail;  // NOLINT
using Catch::Matchers::EndsWith;
using Catch::Matchers::WithinAbs;
using gelex::test::FileFixture;

TEST_CASE("CovarEffectLoader Constructor Tests", "[predictor][covar_effect]")
{
    FileFixture files;

    SECTION(
        "Happy path - valid parameter file with intercept, continuous and "
        "categorical vars")
    {
        auto file_path = files.create_text_file(
            "term\tmean\tstddev\tpercentile_5\tpercentile_95\tess\trhat\n"
            "Intercept\t2.5\t0.1\t2.3\t2.7\t1000\t1.0\n"
            "Age\t0.5\t0.05\t0.4\t0.6\t800\t1.01\n"
            "Height\t-0.2\t0.02\t-0.23\t-0.17\t1200\t1.02\n"
            "Sex_M\t-0.3\t0.02\t-0.33\t-0.27\t1100\t1.01\n"
            "Sex_F\t0.2\t0.02\t0.17\t0.23\t900\t1.03\n"
            "Population_EUR\t0.8\t0.1\t0.6\t1.0\t700\t1.05\n"
            "Population_AFR\t-0.5\t0.08\t-0.63\t-0.37\t850\t1.02\n");

        REQUIRE_NOTHROW(
            [&]()
            {
                CovarEffectLoader loader(file_path);
                const auto& covar_effect = loader.effects();
                REQUIRE(covar_effect.intercept == 2.5);

                const auto& continuous = covar_effect.continuous_coeffs;
                REQUIRE(continuous.size() == 2);
                REQUIRE(continuous.at("Age") == 0.5);
                REQUIRE(continuous.at("Height") == -0.2);

                const auto& categorical = covar_effect.categorical_coeffs;
                REQUIRE(categorical.size() == 2);
                REQUIRE(categorical.at("Sex").size() == 2);
                REQUIRE(categorical.at("Sex").at("M") == -0.3);
                REQUIRE(categorical.at("Sex").at("F") == 0.2);
                REQUIRE(categorical.at("Population").size() == 2);
                REQUIRE(categorical.at("Population").at("EUR") == 0.8);
                REQUIRE(categorical.at("Population").at("AFR") == -0.5);
            }());
    }

    SECTION("Happy path - file with only intercept")
    {
        auto file_path = files.create_text_file(
            "term\tmean\tstddev\tpercentile_5\tpercentile_95\tess\trhat\n"
            "Intercept\t1.8\t0.2\t1.4\t2.2\t1500\t1.01\n");

        REQUIRE_NOTHROW(
            [&]()
            {
                CovarEffectLoader loader(file_path);
                const auto& covar_effect = loader.effects();
                REQUIRE(covar_effect.intercept == 1.8);
                REQUIRE(covar_effect.continuous_coeffs.empty());
                REQUIRE(covar_effect.categorical_coeffs.empty());
            }());
    }

    SECTION("Happy path - file with only continuous variables")
    {
        auto file_path = files.create_text_file(
            "term\tmean\tstddev\tpercentile_5\tpercentile_95\tess\trhat\n"
            "Intercept\t0.0\t0.1\t-0.2\t0.2\t1000\t1.0\n"
            "BMI\t0.3\t0.05\t0.2\t0.4\t800\t1.01\n"
            "Cholesterol\t0.1\t0.03\t0.04\t0.16\t900\t1.02\n");

        REQUIRE_NOTHROW(
            [&]()
            {
                CovarEffectLoader loader(file_path);
                const auto& covar_effect = loader.effects();
                REQUIRE(covar_effect.intercept == 0.0);

                const auto& continuous = covar_effect.continuous_coeffs;
                REQUIRE(continuous.size() == 2);
                REQUIRE(continuous.at("BMI") == 0.3);
                REQUIRE(continuous.at("Cholesterol") == 0.1);
                REQUIRE(covar_effect.categorical_coeffs.empty());
            }());
    }

    SECTION("Happy path - file with only categorical variables")
    {
        auto file_path = files.create_text_file(
            "term\tmean\tstddev\tpercentile_5\tpercentile_95\tess\trhat\n"
            "Intercept\t3.0\t0.3\t2.4\t3.6\t1200\t1.0\n"
            "Group_A\t0.5\t0.1\t0.3\t0.7\t800\t1.01\n"
            "Group_B\t-0.2\t0.08\t-0.34\t-0.06\t700\t1.02\n"
            "Treatment_Placebo\t0.0\t0.05\t-0.1\t0.1\t900\t1.03\n"
            "Treatment_Drug\t0.8\t0.15\t0.5\t1.1\t850\t1.04\n");

        REQUIRE_NOTHROW(
            [&]()
            {
                CovarEffectLoader loader(file_path);
                REQUIRE(loader.effects().intercept == 3.0);
                REQUIRE(loader.effects().continuous_coeffs.empty());

                const auto& categorical = loader.effects().categorical_coeffs;
                REQUIRE(categorical.size() == 2);
                REQUIRE(categorical.at("Group").size() == 2);
                REQUIRE(categorical.at("Group").at("A") == 0.5);
                REQUIRE(categorical.at("Group").at("B") == -0.2);
                REQUIRE(categorical.at("Treatment").size() == 2);
                REQUIRE(categorical.at("Treatment").at("Placebo") == 0.0);
                REQUIRE(categorical.at("Treatment").at("Drug") == 0.8);
            }());
    }

    SECTION("Exception - missing intercept term")
    {
        auto file_path = files.create_text_file(
            "term\tmean\tstddev\tpercentile_5\tpercentile_95\tess\trhat\n"
            "Age\t0.5\t0.05\t0.4\t0.6\t800\t1.01\n"
            "Height\t-0.2\t0.02\t-0.23\t-0.17\t1200\t1.02\n");

        REQUIRE_THROWS_AS(
            CovarEffectLoader(file_path), gelex::FileFormatException);
    }

    SECTION("Exception - empty file")
    {
        auto file_path = files.create_empty_file(".txt");

        REQUIRE_THROWS_AS(
            CovarEffectLoader(file_path), gelex::FileFormatException);
    }

    SECTION("Exception - file with only header line")
    {
        auto file_path = files.create_text_file(
            "term\tmean\tstddev\tpercentile_5\tpercentile_95\tess\trhat\n");

        REQUIRE_THROWS_AS(
            CovarEffectLoader(file_path), gelex::FileFormatException);
    }

    SECTION("Happy path - malformed lines are skipped")
    {
        auto file_path = files.create_text_file(
            "term\tmean\tstddev\tpercentile_5\tpercentile_95\tess\trhat\n"
            "Intercept\t2.0\t0.1\t1.8\t2.2\t1000\t1.0\n"
            "InvalidLine\n"
            "Age\t0.5\t0.05\t0.4\t0.6\t800\t1.01\n"
            "AnotherBadLine\tnot_a_number\n"
            "Sex_M\t-0.3\t0.02\t-0.33\t-0.27\t1100\t1.01\n");

        REQUIRE_NOTHROW(
            [&]()
            {
                CovarEffectLoader loader(file_path);
                const auto& covar_effect = loader.effects();
                REQUIRE(covar_effect.intercept == 2.0);

                const auto& continuous = covar_effect.continuous_coeffs;
                REQUIRE(continuous.size() == 1);
                REQUIRE(continuous.at("Age") == 0.5);

                const auto& categorical = covar_effect.categorical_coeffs;
                REQUIRE(categorical.size() == 1);
                REQUIRE(categorical.at("Sex").size() == 1);
                REQUIRE(categorical.at("Sex").at("M") == -0.3);
            }());
    }

    SECTION("Happy path - duplicate variable names overwrite previous values")
    {
        auto file_path = files.create_text_file(
            "term\tmean\tstddev\tpercentile_5\tpercentile_95\tess\trhat\n"
            "Intercept\t1.0\t0.1\t0.8\t1.2\t1000\t1.0\n"
            "Age\t0.3\t0.05\t0.2\t0.4\t800\t1.01\n"
            "Age\t0.5\t0.05\t0.4\t0.6\t900\t1.02\n"  // Should overwrite
                                                     // previous
            "Group_A\t0.2\t0.03\t0.14\t0.26\t700\t1.03\n"
            "Group_A\t0.3\t0.04\t0.22\t0.38\t750\t1.04\n");  // Should overwrite

        REQUIRE_NOTHROW(
            [&]()
            {
                CovarEffectLoader loader(file_path);

                const auto& covar_effect = loader.effects();

                REQUIRE(covar_effect.intercept == 1.0);

                const auto& continuous = covar_effect.continuous_coeffs;
                REQUIRE(continuous.size() == 1);
                REQUIRE(continuous.at("Age") == 0.5);  // Last value wins

                const auto& categorical = covar_effect.categorical_coeffs;
                REQUIRE(categorical.size() == 1);
                REQUIRE(categorical.at("Group").size() == 1);
                REQUIRE(
                    categorical.at("Group").at("A") == 0.3);  // Last value wins
            }());
    }
}

TEST_CASE("CovarEffectLoader Accessor Tests", "[predictor][covar_effect]")
{
    FileFixture files;

    SECTION("intercept() returns correct value")
    {
        auto file_path = files.create_text_file(
            "term\tmean\tstddev\tpercentile_5\tpercentile_95\tess\trhat\n"
            "Intercept\t3.14159\t0.01\t3.12159\t3.16159\t2000\t1.0\n");

        CovarEffectLoader loader(file_path);
        REQUIRE(loader.effects().intercept == 3.14159);
    }

    SECTION("continuous_coeffs() returns correct mapping")
    {
        auto file_path = files.create_text_file(
            "term\tmean\tstddev\tpercentile_5\tpercentile_95\tess\trhat\n"
            "Intercept\t0.0\t0.1\t-0.2\t0.2\t1000\t1.0\n"
            "Var1\t1.5\t0.2\t1.1\t1.9\t800\t1.01\n"
            "Var2\t-2.0\t0.3\t-2.6\t-1.4\t900\t1.02\n"
            "Var3\t0.8\t0.15\t0.5\t1.1\t700\t1.03\n");

        CovarEffectLoader loader(file_path);
        const auto& coeffs = loader.effects().continuous_coeffs;

        REQUIRE(coeffs.size() == 3);
        REQUIRE(coeffs.at("Var1") == 1.5);
        REQUIRE(coeffs.at("Var2") == -2.0);
        REQUIRE(coeffs.at("Var3") == 0.8);

        // Ensure map is ordered (std::map maintains order)
        auto it = coeffs.begin();
        REQUIRE(it->first == "Var1");
        REQUIRE(it->second == 1.5);
        ++it;
        REQUIRE(it->first == "Var2");
        REQUIRE(it->second == -2.0);
        ++it;
        REQUIRE(it->first == "Var3");
        REQUIRE(it->second == 0.8);
    }

    SECTION("categorical_coeffs() returns correct nested mapping")
    {
        auto file_path = files.create_text_file(
            "term\tmean\tstddev\tpercentile_5\tpercentile_95\tess\trhat\n"
            "Intercept\t1.0\t0.1\t0.8\t1.2\t1000\t1.0\n"
            "Group_A\t0.5\t0.1\t0.3\t0.7\t800\t1.01\n"
            "Group_B\t-0.2\t0.08\t-0.34\t-0.06\t700\t1.02\n"
            "Group_C\t0.3\t0.12\t0.06\t0.54\t900\t1.03\n"
            "Treatment_Placebo\t0.0\t0.05\t-0.1\t0.1\t850\t1.04\n"
            "Treatment_Drug\t0.8\t0.15\t0.5\t1.1\t750\t1.05\n");

        CovarEffectLoader loader(file_path);
        const auto& coeffs = loader.effects().categorical_coeffs;

        REQUIRE(coeffs.size() == 2);
        REQUIRE(coeffs.count("Group") == 1);
        REQUIRE(coeffs.count("Treatment") == 1);

        const auto& group_coeffs = coeffs.at("Group");
        REQUIRE(group_coeffs.size() == 3);
        REQUIRE(group_coeffs.at("A") == 0.5);
        REQUIRE(group_coeffs.at("B") == -0.2);
        REQUIRE(group_coeffs.at("C") == 0.3);

        const auto& treatment_coeffs = coeffs.at("Treatment");
        REQUIRE(treatment_coeffs.size() == 2);
        REQUIRE(treatment_coeffs.at("Placebo") == 0.0);
        REQUIRE(treatment_coeffs.at("Drug") == 0.8);
    }
}

TEST_CASE("CovarEffectLoader Error Handling Tests", "[predictor][covar_effect]")
{
    FileFixture files;

    SECTION("Exception - file does not exist")
    {
        auto non_existent_path = files.generate_random_file_path(".txt");

        REQUIRE_THROWS(CovarEffectLoader(non_existent_path));
    }

    SECTION("Happy path - insufficient columns in data line (skipped)")
    {
        auto file_path = files.create_text_file(
            "term\tmean\tstddev\tpercentile_5\tpercentile_"
            "95\tess\trhat\n"
            "Intercept\t1.0\t0.1\t0.8\t1.2\t1000\t1.0\n"
            "Age\t0.5\t0.05\t0.4\t0.6\t800\n"  // Missing rhat column
            "Height\t0.2\t0.03\t0.14\t0.26\t900\t1."
            "02\tExtraColumn\n");  // Extra column

        REQUIRE_NOTHROW(
            [&]()
            {
                CovarEffectLoader loader(file_path);
                REQUIRE(loader.effects().intercept == 1.0);

                // Age line should be skipped due to missing column
                // Height line should parse successfully (extra column ignored)
                const auto& continuous = loader.effects().continuous_coeffs;
                REQUIRE(continuous.size() == 1);
                REQUIRE(continuous.at("Height") == 0.2);
            }());
    }

    SECTION("Happy path - non-numeric mean value (skipped)")
    {
        auto file_path = files.create_text_file(
            "term\tmean\tstddev\tpercentile_5\tpercentile_95\tess\trhat\n"
            "Intercept\t1.0\t0.1\t0.8\t1.2\t1000\t1.0\n"
            "Age\tnot_a_number\t0.05\t0.4\t0.6\t800\t1.01\n"
            "Height\t0.2\t0.03\t0.14\t0.26\t900\t1.02\n");

        REQUIRE_NOTHROW(
            [&]()
            {
                CovarEffectLoader loader(file_path);
                REQUIRE(loader.effects().intercept == 1.0);

                // Age line should be skipped due to non-numeric mean
                const auto& continuous = loader.effects().continuous_coeffs;
                REQUIRE(continuous.size() == 1);
                REQUIRE(continuous.at("Height") == 0.2);
            }());
    }

    SECTION("Happy path - empty values in data line (skipped)")
    {
        auto file_path = files.create_text_file(
            "term\tmean\tstddev\tpercentile_5\tpercentile_95\tess\trhat\n"
            "Intercept\t1.0\t0.1\t0.8\t1.2\t1000\t1.0\n"
            "Age\t\t0.05\t0.4\t0.6\t800\t1.01\n"  // Empty mean
            "Height\t0.2\t0.03\t0.14\t0.26\t900\t1.02\n");

        REQUIRE_NOTHROW(
            [&]()
            {
                CovarEffectLoader loader(file_path);
                REQUIRE(loader.effects().intercept == 1.0);

                // Age line should be skipped due to empty mean
                const auto& continuous = loader.effects().continuous_coeffs;
                REQUIRE(continuous.size() == 1);
                REQUIRE(continuous.at("Height") == 0.2);
            }());
    }
}
