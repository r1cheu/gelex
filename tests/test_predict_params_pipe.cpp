#include <string_view>

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_exception.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>

#include "../src/predictor/predict_params_pipe.h"
#include "file_fixture.h"
#include "gelex/exception.h"

namespace fs = std::filesystem;

using Catch::Matchers::EndsWith;
using gelex::PredictParamsPipe;
using gelex::SnpEffects;
using gelex::detail::CovarEffects;
using gelex::test::FileFixture;

TEST_CASE(
    "PredictParamsPipe - Constructor Success Scenarios",
    "[predictor][predict_params]")
{
    FileFixture files;

    SECTION("Happy path - both valid files loaded successfully")
    {
        std::string snp_content
            = "ID\tA1\tA2\tA1Frq\tAdd\tDom\n"
              "rs001\tA\tC\t0.25\t0.123\t0.045\n"
              "rs002\tT\tG\t0.75\t-0.456\t0.089\n"
              "rs003\tC\tA\t0.50\t0.789\t-0.012\n";

        auto snp_path = files.create_text_file(snp_content, ".snp.eff");

        std::string covar_content
            = "term\tmean\tstddev\tpercentile_5\tpercentile_95\tess\trhat\n"
              "Intercept\t2.5\t0.1\t2.3\t2.7\t1000\t1.0\n"
              "Age\t0.5\t0.05\t0.4\t0.6\t800\t1.01\n"
              "Height\t-0.2\t0.02\t-0.23\t-0.17\t1200\t1.02\n";

        auto covar_path = files.create_text_file(covar_content, ".covar.eff");

        PredictParamsPipe::Config config;
        config.snp_effect_path = snp_path;
        config.covar_effect_path = covar_path;

        PredictParamsPipe pipe(config);

        const auto& snp_effects = pipe.snp_effects();
        REQUIRE(snp_effects.size() == 3);

        const auto& covar_effects = pipe.covar_effects();
        REQUIRE(covar_effects.intercept == 2.5);
        REQUIRE(covar_effects.continuous_coeffs.size() == 2);
        REQUIRE(covar_effects.categorical_coeffs.empty());
    }

    SECTION("Exception - covariate effect file has invalid format (propagated)")
    {
        std::string snp_content
            = "ID\tA1\tA2\tA1Frq\tAdd\tDom\n"
              "rs001\tA\tC\t0.25\t0.123\t0.045\n";

        auto snp_path = files.create_text_file(snp_content, ".snp.eff");

        std::string invalid_covar_content
            = "term\tmean\tstddev\n"  // Missing required columns
              "Intercept\t1.0\t0.1\n";

        auto covar_path
            = files.create_text_file(invalid_covar_content, ".covar.eff");

        PredictParamsPipe::Config config;
        config.snp_effect_path = snp_path;
        config.covar_effect_path = covar_path;

        REQUIRE_THROWS_AS(
            PredictParamsPipe(config), gelex::FileFormatException);
    }
}

TEST_CASE("PredictParamsPipe - Accessor Methods", "[predictor][predict_params]")
{
    FileFixture files;

    SECTION("snp_effects() returns const reference")
    {
        std::string snp_content
            = "ID\tA1\tA2\tA1Frq\tAdd\tDom\n"
              "rs001\tA\tC\t0.25\t0.123\t0.045\n";

        auto snp_path = files.create_text_file(snp_content, ".snp.eff");

        std::string covar_content
            = "term\tmean\tstddev\tpercentile_5\tpercentile_95\tess\trhat\n"
              "Intercept\t2.0\t0.1\t1.8\t2.2\t1000\t1.0\n";

        auto covar_path = files.create_text_file(covar_content, ".covar.eff");

        PredictParamsPipe::Config config;
        config.snp_effect_path = snp_path;
        config.covar_effect_path = covar_path;

        PredictParamsPipe pipe(config);

        const auto& effects1 = pipe.snp_effects();
        const auto& effects2 = pipe.snp_effects();
        REQUIRE(&effects1 == &effects2);  // Same reference

        REQUIRE(effects1.size() == 1);
        REQUIRE(effects1.contains("rs001"));
    }

    SECTION("covar_effects() returns const reference")
    {
        std::string snp_content
            = "ID\tA1\tA2\tA1Frq\tAdd\tDom\n"
              "rs001\tA\tC\t0.25\t0.123\t0.045\n";

        auto snp_path = files.create_text_file(snp_content, ".snp.eff");

        std::string covar_content
            = "term\tmean\tstddev\tpercentile_5\tpercentile_95\tess\trhat\n"
              "Intercept\t3.0\t0.2\t2.6\t3.4\t800\t1.01\n"
              "Age\t0.5\t0.05\t0.4\t0.6\t900\t1.02\n";

        auto covar_path = files.create_text_file(covar_content, ".covar.eff");

        PredictParamsPipe::Config config;
        config.snp_effect_path = snp_path;
        config.covar_effect_path = covar_path;

        PredictParamsPipe pipe(config);

        const auto& effects1 = pipe.covar_effects();
        const auto& effects2 = pipe.covar_effects();
        REQUIRE(&effects1 == &effects2);  // Same reference

        REQUIRE(effects1.intercept == 3.0);
        REQUIRE(effects1.continuous_coeffs.size() == 1);
        REQUIRE(effects1.continuous_coeffs.at("Age") == 0.5);
    }
}

TEST_CASE("PredictParamsPipe - Move Semantics", "[predictor][predict_params]")
{
    FileFixture files;

    SECTION("take_snp_effects() moves data out")
    {
        std::string snp_content
            = "ID\tA1\tA2\tA1Frq\tAdd\tDom\n"
              "rs001\tA\tC\t0.25\t0.123\t0.045\n"
              "rs002\tT\tG\t0.75\t-0.456\t0.089\n";

        auto snp_path = files.create_text_file(snp_content, ".snp.eff");

        std::string covar_content
            = "term\tmean\tstddev\tpercentile_5\tpercentile_95\tess\trhat\n"
              "Intercept\t2.0\t0.1\t1.8\t2.2\t1000\t1.0\n";

        auto covar_path = files.create_text_file(covar_content, ".covar.eff");

        PredictParamsPipe::Config config;
        config.snp_effect_path = snp_path;
        config.covar_effect_path = covar_path;

        PredictParamsPipe pipe(config);
        REQUIRE(pipe.snp_effects().size() == 2);

        SnpEffects moved_effects = std::move(pipe).take_snp_effects();
        REQUIRE(moved_effects.size() == 2);
        REQUIRE(pipe.snp_effects().empty());  // Data moved out
    }

    SECTION("take_covar_effects() moves data out")
    {
        std::string snp_content
            = "ID\tA1\tA2\tA1Frq\tAdd\tDom\n"
              "rs001\tA\tC\t0.25\t0.123\t0.045\n";

        auto snp_path = files.create_text_file(snp_content, ".snp.eff");

        std::string covar_content
            = "term\tmean\tstddev\tpercentile_5\tpercentile_95\tess\trhat\n"
              "Intercept\t3.0\t0.2\t2.6\t3.4\t800\t1.01\n"
              "Age\t0.5\t0.05\t0.4\t0.6\t900\t1.02\n";

        auto covar_path = files.create_text_file(covar_content, ".covar.eff");

        PredictParamsPipe::Config config;
        config.snp_effect_path = snp_path;
        config.covar_effect_path = covar_path;

        PredictParamsPipe pipe(config);
        REQUIRE(pipe.covar_effects().intercept == 3.0);

        CovarEffects moved_effects = std::move(pipe).take_covar_effects();
        REQUIRE(moved_effects.intercept == 3.0);
        REQUIRE(moved_effects.continuous_coeffs.size() == 1);

        // After move, covar_effects maps are empty (data moved out)
        const auto& after_effects = pipe.covar_effects();
        REQUIRE(after_effects.continuous_coeffs.empty());
        REQUIRE(after_effects.categorical_coeffs.empty());
    }
}

TEST_CASE("PredictParamsPipe - Edge Cases", "[predictor][predict_params]")
{
    FileFixture files;

    SECTION("Empty SNP effect file (only header)")
    {
        std::string snp_content = "ID\tA1\tA2\tA1Frq\tAdd\tDom\n";
        auto snp_path = files.create_text_file(snp_content, ".snp.eff");

        std::string covar_content
            = "term\tmean\tstddev\tpercentile_5\tpercentile_95\tess\trhat\n"
              "Intercept\t1.0\t0.1\t0.8\t1.2\t1000\t1.0\n";

        auto covar_path = files.create_text_file(covar_content, ".covar.eff");

        PredictParamsPipe::Config config;
        config.snp_effect_path = snp_path;
        config.covar_effect_path = covar_path;

        PredictParamsPipe pipe(config);

        REQUIRE(pipe.snp_effects().empty());
        REQUIRE(pipe.covar_effects().intercept == 1.0);
    }

    SECTION("Covariate file with only intercept")
    {
        std::string snp_content
            = "ID\tA1\tA2\tA1Frq\tAdd\tDom\n"
              "rs001\tA\tC\t0.25\t0.123\t0.045\n";

        auto snp_path = files.create_text_file(snp_content, ".snp.eff");

        std::string covar_content
            = "term\tmean\tstddev\tpercentile_5\tpercentile_95\tess\trhat\n"
              "Intercept\t0.0\t0.1\t-0.2\t0.2\t1000\t1.0\n";

        auto covar_path = files.create_text_file(covar_content, ".covar.eff");

        PredictParamsPipe::Config config;
        config.snp_effect_path = snp_path;
        config.covar_effect_path = covar_path;

        PredictParamsPipe pipe(config);

        REQUIRE(pipe.snp_effects().size() == 1);
        REQUIRE(pipe.covar_effects().intercept == 0.0);
        REQUIRE(pipe.covar_effects().continuous_coeffs.empty());
        REQUIRE(pipe.covar_effects().categorical_coeffs.empty());
    }
}

TEST_CASE(
    "PredictParamsPipe - Constructor Error Scenarios",
    "[predictor][predict_params]")
{
    FileFixture files;

    SECTION("Exception - empty SNP effect path with message")
    {
        std::string covar_content
            = "term\tmean\tstddev\tpercentile_5\tpercentile_95\tess\trhat\n"
              "Intercept\t1.0\t0.1\t0.8\t1.2\t1000\t1.0\n";

        auto covar_path = files.create_text_file(covar_content, ".covar.eff");

        PredictParamsPipe::Config config;
        config.snp_effect_path = "";
        config.covar_effect_path = covar_path;

        REQUIRE_THROWS_AS(
            PredictParamsPipe(config), gelex::InvalidInputException);
    }

    SECTION("Exception - empty covariate effect path with message")
    {
        std::string snp_content
            = "ID\tA1\tA2\tA1Frq\tAdd\tDom\n"
              "rs001\tA\tC\t0.25\t0.123\t0.045\n";

        auto snp_path = files.create_text_file(snp_content, ".snp.eff");

        PredictParamsPipe::Config config;
        config.snp_effect_path = snp_path;
        config.covar_effect_path = "";

        REQUIRE_THROWS_AS(
            PredictParamsPipe(config), gelex::InvalidInputException);
    }
}
