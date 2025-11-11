
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <filesystem>
#include <fstream>

#include "../src/predictor/covariate_processor.h"

using namespace gelex::detail;
using namespace Catch::Matchers;

class CovariateProcessorTestFixture
{
   public:
    CovariateProcessorTestFixture()
    {
        // Create a temporary parameter file for testing
        test_param_file_
            = std::filesystem::temp_directory_path() / "test_param.param";

        std::ofstream file(test_param_file_);
        file << "term\tmean\tstddev\t5%\t95%\tess\trhat\n";
        file << "Intercept\t81.7907\t0.692097\t80.6663\t82.8195\t29.807\t1."
                "04149\n";
        file << "Group_A\t8.1234\t0.512345\t7.2345\t9.0123\t45.678\t1.00234\n";
        file << "Group_B\t-5.6789\t0.456789\t-6.5432\t-4.5678\t39.012\t1."
                "00012\n";
        file << "Age\t0.5\t0.1\t0.3\t0.7\t50.0\t1.001\n";
        file << "Height\t2.3\t0.2\t2.0\t2.6\t40.0\t1.002\n";
        file.close();
    }

    ~CovariateProcessorTestFixture()
    {
        // Clean up the temporary file
        if (std::filesystem::exists(test_param_file_))
        {
            std::filesystem::remove(test_param_file_);
        }
    }

    std::filesystem::path test_param_file_;
};

TEST_CASE_METHOD(
    CovariateProcessorTestFixture,
    "CovariateProcessor handles continuous variables",
    "[covariate_processor]")
{
    auto processor = CovariateProcessor::create(test_param_file_);
    REQUIRE(processor.has_value());

    IndividualData data;
    data.continuous_values["Age"] = 30.0;
    data.continuous_values["Height"] = 170.0;

    double prediction = processor->predict(data);
    double expected = 81.7907 + (30.0 * 0.5) + (170.0 * 2.3);

    REQUIRE_THAT(prediction, WithinAbs(expected, 1e-6));
}

TEST_CASE_METHOD(
    CovariateProcessorTestFixture,
    "CovariateProcessor handles categorical variables",
    "[covariate_processor]")
{
    auto processor = CovariateProcessor::create(test_param_file_);
    REQUIRE(processor.has_value());

    IndividualData data;
    data.categorical_values["Group"] = "A";

    double prediction = processor->predict(data);
    double expected = 81.7907 + 8.1234;

    REQUIRE_THAT(prediction, WithinAbs(expected, 1e-6));
}

TEST_CASE_METHOD(
    CovariateProcessorTestFixture,
    "CovariateProcessor handles mixed variables",
    "[covariate_processor]")
{
    auto processor = CovariateProcessor::create(test_param_file_);
    REQUIRE(processor.has_value());

    IndividualData data;
    data.continuous_values["Age"] = 25.0;
    data.categorical_values["Group"] = "B";

    double prediction = processor->predict(data);
    double expected = 81.7907 + (25.0 * 0.5) + (-5.6789);

    REQUIRE_THAT(prediction, WithinAbs(expected, 1e-6));
}

TEST_CASE_METHOD(
    CovariateProcessorTestFixture,
    "CovariateProcessor ignores unknown variables",
    "[covariate_processor]")
{
    auto processor = CovariateProcessor::create(test_param_file_);
    REQUIRE(processor.has_value());

    IndividualData data;
    data.continuous_values["UnknownVar"] = 100.0;
    data.categorical_values["UnknownCat"] = "SomeValue";

    double prediction = processor->predict(data);
    double expected = 81.7907;  // Only intercept

    REQUIRE_THAT(prediction, WithinAbs(expected, 1e-6));
}

TEST_CASE_METHOD(
    CovariateProcessorTestFixture,
    "CovariateProcessor handles empty data",
    "[covariate_processor]")
{
    auto processor = CovariateProcessor::create(test_param_file_);
    REQUIRE(processor.has_value());

    IndividualData data;  // Empty data

    double prediction = processor->predict(data);
    double expected = 81.7907;  // Only intercept

    REQUIRE_THAT(prediction, WithinAbs(expected, 1e-6));
}

TEST_CASE(
    "CovariateProcessor returns error for non-existent file",
    "[covariate_processor]")
{
    auto processor = CovariateProcessor::create("non_existent_file.param");
    REQUIRE(!processor.has_value());
    REQUIRE(processor.error().code == gelex::ErrorCode::FileNotFound);
}

TEST_CASE(
    "CovariateProcessor handles malformed parameter file",
    "[covariate_processor]")
{
    // Create a temporary malformed parameter file
    std::filesystem::path malformed_file
        = std::filesystem::temp_directory_path() / "malformed.param";

    {
        std::ofstream file(malformed_file);
        file << "term\tmean\tstddev\t5%\t95%\tess\trhat\n";
        file << "Intercept\t81.7907\t0.692097\t80.6663\t82.8195\t29.807\t1."
                "04149\n";
        file << "InvalidLine\n";  // Malformed line
        file << "Group_A\t8.1234\t0.512345\t7.2345\t9.0123\t45.678\t1.00234\n";
        file.close();
    }

    // Should succeed, should skip malformed lines
    auto processor = CovariateProcessor::create(malformed_file);
    REQUIRE(processor.has_value());

    // Clean up
    std::filesystem::remove(malformed_file);
}

TEST_CASE(
    "CovariateProcessor returns error when no intercept found",
    "[covariate_processor]")
{
    // Create a parameter file without intercept
    std::filesystem::path no_intercept_file
        = std::filesystem::temp_directory_path() / "no_intercept.param";

    {
        std::ofstream file(no_intercept_file);
        file << "term\tmean\tstddev\t5%\t95%\tess\trhat\n";
        file << "Group_A\t8.1234\t0.512345\t7.2345\t9.0123\t45.678\t1.00234\n";
        file << "Age\t0.5\t0.1\t0.3\t0.7\t50.0\t1.001\n";
        file.close();
    }

    auto processor = CovariateProcessor::create(no_intercept_file);
    REQUIRE(!processor.has_value());
    REQUIRE(processor.error().code == gelex::ErrorCode::InvalidData);

    // Clean up
    std::filesystem::remove(no_intercept_file);
}
