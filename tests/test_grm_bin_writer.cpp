#include <cmath>
#include <cstddef>
#include <filesystem>
#include <fstream>
#include <limits>
#include <vector>

#include <Eigen/Core>
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_exception.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>

#include "../src/data/grm_bin_writer.h"
#include "file_fixture.h"
#include "gelex/exception.h"

namespace fs = std::filesystem;

using namespace gelex::detail;  // NOLINT
using Catch::Matchers::ContainsSubstring;
using Catch::Matchers::MessageMatches;
using gelex::test::FileFixture;

// Helper function to calculate expected file size for lower triangle
auto expected_file_size(Eigen::Index n) -> size_t
{
    return static_cast<size_t>(n * (n + 1) / 2) * sizeof(float);
}

// Helper function to read lower triangle values from binary file
auto read_lower_triangle(const fs::path& file_path, Eigen::Index n)
    -> std::vector<float>
{
    std::ifstream file(file_path, std::ios::binary);
    const auto num_elements = static_cast<size_t>(n * (n + 1) / 2);
    std::vector<float> values(num_elements);
    file.read(
        reinterpret_cast<char*>(values.data()),
        static_cast<std::streamsize>(num_elements * sizeof(float)));
    return values;
}

// ============================================================================
// Constructor tests
// ============================================================================

TEST_CASE(
    "GrmBinWriter - Constructor and path access",
    "[grm_bin_writer][construction]")
{
    FileFixture files;

    SECTION("Happy path - create writer with valid path")
    {
        auto file_path = files.generate_random_file_path(".grm.bin");

        REQUIRE_NOTHROW(
            [&]()
            {
                GrmBinWriter writer(file_path);
                REQUIRE(writer.path() == file_path);
            }());
    }
}

// ============================================================================
// Empty matrix tests
// ============================================================================

TEST_CASE("GrmBinWriter - Write empty matrix", "[grm_bin_writer][empty]")
{
    FileFixture files;

    SECTION("Happy path - write empty matrix (0x0) produces empty file")
    {
        auto file_path = files.generate_random_file_path(".grm.bin");

        Eigen::MatrixXd empty_matrix(0, 0);

        REQUIRE_NOTHROW(
            [&]()
            {
                GrmBinWriter writer(file_path);
                writer.write(empty_matrix);
            }());

        REQUIRE(fs::exists(file_path));
        REQUIRE(fs::file_size(file_path) == 0);
    }
}

// ============================================================================
// Basic write tests
// ============================================================================

TEST_CASE("GrmBinWriter - Write 1x1 matrix", "[grm_bin_writer][basic]")
{
    FileFixture files;

    SECTION("Happy path - write single element matrix")
    {
        auto file_path = files.generate_random_file_path(".grm.bin");

        Eigen::MatrixXd matrix(1, 1);
        matrix << 1.5;

        REQUIRE_NOTHROW(
            [&]()
            {
                GrmBinWriter writer(file_path);
                writer.write(matrix);
            }());

        REQUIRE(fs::exists(file_path));
        REQUIRE(fs::file_size(file_path) == expected_file_size(1));

        auto values = read_lower_triangle(file_path, 1);
        REQUIRE(values.size() == 1);
        REQUIRE(values[0] == 1.5F);
    }
}

TEST_CASE("GrmBinWriter - Write 3x3 matrix", "[grm_bin_writer][basic]")
{
    FileFixture files;

    SECTION("Happy path - verify lower triangle order")
    {
        auto file_path = files.generate_random_file_path(".grm.bin");

        // Create a 3x3 matrix with distinct values in lower triangle
        // clang-format off
        Eigen::MatrixXd matrix(3, 3);
        matrix << 1.0, 0.0, 0.0,
                  2.0, 3.0, 0.0,
                  4.0, 5.0, 6.0;
        // clang-format on

        REQUIRE_NOTHROW(
            [&]()
            {
                GrmBinWriter writer(file_path);
                writer.write(matrix);
            }());

        REQUIRE(fs::exists(file_path));
        // Lower triangle has 3*(3+1)/2 = 6 elements
        REQUIRE(fs::file_size(file_path) == expected_file_size(3));

        auto values = read_lower_triangle(file_path, 3);
        REQUIRE(values.size() == 6);

        // Expected order: (0,0), (1,0), (1,1), (2,0), (2,1), (2,2)
        // Values:           1.0,   2.0,   3.0,   4.0,   5.0,   6.0
        REQUIRE(values[0] == 1.0F);  // (0,0)
        REQUIRE(values[1] == 2.0F);  // (1,0)
        REQUIRE(values[2] == 3.0F);  // (1,1)
        REQUIRE(values[3] == 4.0F);  // (2,0)
        REQUIRE(values[4] == 5.0F);  // (2,1)
        REQUIRE(values[5] == 6.0F);  // (2,2)
    }
}

TEST_CASE("GrmBinWriter - Write medium matrix", "[grm_bin_writer][basic]")
{
    FileFixture files;

    SECTION("Happy path - write 10x10 matrix")
    {
        auto file_path = files.generate_random_file_path(".grm.bin");

        const Eigen::Index n = 10;
        Eigen::MatrixXd matrix = Eigen::MatrixXd::Random(n, n);
        // Make it symmetric for realistic GRM
        matrix = (matrix + matrix.transpose()) / 2.0;

        REQUIRE_NOTHROW(
            [&]()
            {
                GrmBinWriter writer(file_path);
                writer.write(matrix);
            }());

        REQUIRE(fs::exists(file_path));
        REQUIRE(fs::file_size(file_path) == expected_file_size(n));

        auto values = read_lower_triangle(file_path, n);
        REQUIRE(values.size() == static_cast<size_t>(n * (n + 1) / 2));

        // Verify all lower triangle elements in correct order
        size_t idx = 0;
        for (Eigen::Index i = 0; i < n; ++i)
        {
            for (Eigen::Index j = 0; j <= i; ++j)
            {
                REQUIRE(values[idx] == static_cast<float>(matrix(i, j)));
                ++idx;
            }
        }
    }
}

// ============================================================================
// Numerical verification tests
// ============================================================================

TEST_CASE(
    "GrmBinWriter - Lower triangle order verification",
    "[grm_bin_writer][numerical]")
{
    FileFixture files;

    SECTION("Happy path - verify exact element ordering for 4x4 matrix")
    {
        auto file_path = files.generate_random_file_path(".grm.bin");

        // Create a 4x4 matrix where each element value encodes its position
        // value = i * 10 + j (e.g., element (2,1) has value 21)
        Eigen::MatrixXd matrix(4, 4);
        for (Eigen::Index i = 0; i < 4; ++i)
        {
            for (Eigen::Index j = 0; j < 4; ++j)
            {
                matrix(i, j) = static_cast<double>(i * 10 + j);
            }
        }

        REQUIRE_NOTHROW(
            [&]()
            {
                GrmBinWriter writer(file_path);
                writer.write(matrix);
            }());

        auto values = read_lower_triangle(file_path, 4);

        // Expected order: (0,0), (1,0), (1,1), (2,0), (2,1), (2,2), (3,0),
        // (3,1), (3,2), (3,3) Values:     0,    10,    11,    20,    21, 22,
        // 30,    31,    32,    33
        REQUIRE(values.size() == 10);
        REQUIRE(values[0] == 0.0F);   // (0,0)
        REQUIRE(values[1] == 10.0F);  // (1,0)
        REQUIRE(values[2] == 11.0F);  // (1,1)
        REQUIRE(values[3] == 20.0F);  // (2,0)
        REQUIRE(values[4] == 21.0F);  // (2,1)
        REQUIRE(values[5] == 22.0F);  // (2,2)
        REQUIRE(values[6] == 30.0F);  // (3,0)
        REQUIRE(values[7] == 31.0F);  // (3,1)
        REQUIRE(values[8] == 32.0F);  // (3,2)
        REQUIRE(values[9] == 33.0F);  // (3,3)
    }
}

TEST_CASE(
    "GrmBinWriter - Double to float conversion",
    "[grm_bin_writer][numerical]")
{
    FileFixture files;

    SECTION("Happy path - verify precision reduction from double to float")
    {
        auto file_path = files.generate_random_file_path(".grm.bin");

        // Use values that demonstrate float precision limits
        // float has ~7 significant decimal digits
        Eigen::MatrixXd matrix(2, 2);
        matrix << 1.23456789012345, 0.0, 9.87654321098765, 0.00000012345678;

        REQUIRE_NOTHROW(
            [&]()
            {
                GrmBinWriter writer(file_path);
                writer.write(matrix);
            }());

        auto values = read_lower_triangle(file_path, 2);
        REQUIRE(values.size() == 3);

        // Verify values are converted to float (reduced precision)
        REQUIRE(values[0] == static_cast<float>(1.23456789012345));
        REQUIRE(values[1] == static_cast<float>(9.87654321098765));
        REQUIRE(values[2] == static_cast<float>(0.00000012345678));

        // Verify that float precision is less than double
        // The original double values should not exactly match after conversion
        constexpr double d0 = 1.23456789012345;
        constexpr double d1 = 9.87654321098765;
        REQUIRE(static_cast<double>(values[0]) != d0);
        REQUIRE(static_cast<double>(values[1]) != d1);
    }
}

// ============================================================================
// Exception tests
// ============================================================================

TEST_CASE(
    "GrmBinWriter - Non-square matrix exception",
    "[grm_bin_writer][exception]")
{
    FileFixture files;

    SECTION("Exception - non-square matrix (3x5) throws InvalidInputException")
    {
        auto file_path = files.generate_random_file_path(".grm.bin");

        Eigen::MatrixXd matrix(3, 5);
        matrix.setOnes();

        GrmBinWriter writer(file_path);
        REQUIRE_THROWS_AS(writer.write(matrix), gelex::InvalidInputException);
    }

    SECTION("Exception - non-square matrix (5x3) throws InvalidInputException")
    {
        auto file_path = files.generate_random_file_path(".grm.bin");

        Eigen::MatrixXd matrix(5, 3);
        matrix.setOnes();

        GrmBinWriter writer(file_path);
        REQUIRE_THROWS_AS(writer.write(matrix), gelex::InvalidInputException);
    }

    SECTION("Exception - error message contains dimensions")
    {
        auto file_path = files.generate_random_file_path(".grm.bin");

        Eigen::MatrixXd matrix(3, 7);
        matrix.setOnes();

        GrmBinWriter writer(file_path);
        REQUIRE_THROWS_MATCHES(
            writer.write(matrix),
            gelex::InvalidInputException,
            MessageMatches(ContainsSubstring("3x7")));
    }
}

// ============================================================================
// Special values tests
// ============================================================================

TEST_CASE(
    "GrmBinWriter - Write matrix with special values",
    "[grm_bin_writer][special]")
{
    FileFixture files;

    SECTION("Happy path - write matrix with inf values")
    {
        auto file_path = files.generate_random_file_path(".grm.bin");

        const auto inf = std::numeric_limits<double>::infinity();
        const auto neg_inf = -std::numeric_limits<double>::infinity();

        Eigen::MatrixXd matrix(2, 2);
        matrix << inf, 0.0, 1.0, neg_inf;

        REQUIRE_NOTHROW(
            [&]()
            {
                GrmBinWriter writer(file_path);
                writer.write(matrix);
            }());

        auto values = read_lower_triangle(file_path, 2);
        REQUIRE(values.size() == 3);

        REQUIRE(std::isinf(values[0]));
        REQUIRE(values[0] > 0);  // +inf
        REQUIRE(values[1] == 1.0F);
        REQUIRE(std::isinf(values[2]));
        REQUIRE(values[2] < 0);  // -inf
    }

    SECTION("Happy path - write matrix with NaN values")
    {
        auto file_path = files.generate_random_file_path(".grm.bin");

        const auto nan = std::numeric_limits<double>::quiet_NaN();

        Eigen::MatrixXd matrix(2, 2);
        matrix << nan, 0.0, 2.5, nan;

        REQUIRE_NOTHROW(
            [&]()
            {
                GrmBinWriter writer(file_path);
                writer.write(matrix);
            }());

        auto values = read_lower_triangle(file_path, 2);
        REQUIRE(values.size() == 3);

        REQUIRE(std::isnan(values[0]));
        REQUIRE(values[1] == 2.5F);
        REQUIRE(std::isnan(values[2]));
    }

    SECTION("Happy path - write matrix with very small values")
    {
        auto file_path = files.generate_random_file_path(".grm.bin");

        // Use values near float's minimum positive value
        const auto small_val
            = static_cast<double>(std::numeric_limits<float>::min());

        Eigen::MatrixXd matrix(2, 2);
        matrix << small_val, 0.0, small_val / 2.0, small_val * 2.0;

        REQUIRE_NOTHROW(
            [&]()
            {
                GrmBinWriter writer(file_path);
                writer.write(matrix);
            }());

        auto values = read_lower_triangle(file_path, 2);
        REQUIRE(values.size() == 3);

        REQUIRE(values[0] == static_cast<float>(small_val));
        // small_val / 2.0 may become subnormal or zero
        REQUIRE(values[2] == static_cast<float>(small_val * 2.0));
    }
}

// ============================================================================
// Buffer size verification
// ============================================================================

TEST_CASE("GrmBinWriter - Buffer size verification", "[grm_bin_writer][buffer]")
{
    FileFixture files;

    SECTION("Happy path - verify default buffer size constant")
    {
        REQUIRE(
            GrmBinWriter::kDefaultBufferSize == static_cast<size_t>(64 * 1024));
    }

    SECTION("Happy path - write with default buffer should work")
    {
        auto file_path = files.generate_random_file_path(".grm.bin");

        // Create a reasonably sized matrix
        const Eigen::Index n = 50;
        Eigen::MatrixXd matrix = Eigen::MatrixXd::Random(n, n);
        matrix = (matrix + matrix.transpose()) / 2.0;

        REQUIRE_NOTHROW(
            [&]()
            {
                GrmBinWriter writer(file_path);
                writer.write(matrix);
            }());

        REQUIRE(fs::exists(file_path));
        REQUIRE(fs::file_size(file_path) == expected_file_size(n));
    }
}

// ============================================================================
// File size verification tests
// ============================================================================

TEST_CASE(
    "GrmBinWriter - File size verification",
    "[grm_bin_writer][file_size]")
{
    FileFixture files;

    SECTION("Happy path - file size for n=5")
    {
        auto file_path = files.generate_random_file_path(".grm.bin");
        const Eigen::Index n = 5;
        Eigen::MatrixXd matrix = Eigen::MatrixXd::Random(n, n);

        {
            GrmBinWriter writer(file_path);
            writer.write(matrix);
        }

        REQUIRE(fs::file_size(file_path) == expected_file_size(n));
    }

    SECTION("Happy path - file size for n=10")
    {
        auto file_path = files.generate_random_file_path(".grm.bin");
        const Eigen::Index n = 10;
        Eigen::MatrixXd matrix = Eigen::MatrixXd::Random(n, n);

        {
            GrmBinWriter writer(file_path);
            writer.write(matrix);
        }

        REQUIRE(fs::file_size(file_path) == expected_file_size(n));
    }

    SECTION("Happy path - file size for n=20")
    {
        auto file_path = files.generate_random_file_path(".grm.bin");
        const Eigen::Index n = 20;
        Eigen::MatrixXd matrix = Eigen::MatrixXd::Random(n, n);

        {
            GrmBinWriter writer(file_path);
            writer.write(matrix);
        }

        REQUIRE(fs::file_size(file_path) == expected_file_size(n));
    }
}
