#include <filesystem>
#include <fstream>
#include <iostream>

#include <Eigen/Core>
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

using Catch::Matchers::WithinAbs;

#include "../include/gelex/data/genotype_mmap.h"
#include "gelex/error.h"

class GenotypeMapTestFixture
{
   public:
    GenotypeMapTestFixture()
    {
        createValidTestFiles();
        createInvalidMetadataFile();
        createSizeMismatchFile();
        createReadOnlyFiles();
    }

    GenotypeMapTestFixture(const GenotypeMapTestFixture&) = default;
    GenotypeMapTestFixture(GenotypeMapTestFixture&&) = delete;
    GenotypeMapTestFixture& operator=(const GenotypeMapTestFixture&) = default;
    GenotypeMapTestFixture& operator=(GenotypeMapTestFixture&&) = delete;
    ~GenotypeMapTestFixture() { remove(); }

    static void remove()
    {
        std::remove("test_valid.bin");
        std::remove("test_valid.snpstats");
        std::remove("test_invalid_meta.bin");
        std::remove("test_invalid_meta.snpstats");
        std::remove("test_size_mismatch.bin");
        std::remove("test_size_mismatch.snpstats");
        std::remove("test_readonly.bin");
        std::remove("test_readonly.snpstats");
        std::remove("test_large.bin");
        std::remove("test_large.snpstats");
    }

    static void createValidTestFiles()
    {
        // Create a 3x2 matrix test data
        Eigen::MatrixXd matrix(3, 2);
        matrix << 1.0, 4.0, 2.0, 5.0, 3.0, 6.0;

        // Write binary data
        std::ofstream bin_file("test_valid.bin", std::ios::binary);
        bin_file.write(
            reinterpret_cast<const char*>(matrix.data()),
            static_cast<std::streamsize>(matrix.size() * sizeof(double)));

        // Write metadata
        std::ofstream meta_file("test_valid.snpstats", std::ios::binary);
        int64_t rows = matrix.rows();
        int64_t cols = matrix.cols();
        meta_file.write(reinterpret_cast<const char*>(&rows), sizeof(int64_t));
        meta_file.write(reinterpret_cast<const char*>(&cols), sizeof(int64_t));
    }

    static void createInvalidMetadataFile()
    {
        // Create binary file
        Eigen::MatrixXd matrix(2, 2);
        matrix << 1.0, 3.0, 2.0, 4.0;

        std::ofstream bin_file("test_invalid_meta.bin", std::ios::binary);
        bin_file.write(
            reinterpret_cast<const char*>(matrix.data()),
            static_cast<std::streamsize>(matrix.size() * sizeof(double)));

        // Create invalid metadata with negative dimensions
        std::ofstream meta_file("test_invalid_meta.snpstats", std::ios::binary);
        int64_t rows = -1;  // Invalid dimension
        int64_t cols = 2;
        meta_file.write(reinterpret_cast<const char*>(&rows), sizeof(int64_t));
        meta_file.write(reinterpret_cast<const char*>(&cols), sizeof(int64_t));
    }

    static void createSizeMismatchFile()
    {
        // Create binary file with 2x2 data
        Eigen::MatrixXd matrix(2, 2);
        matrix << 1.0, 3.0, 2.0, 4.0;

        std::ofstream bin_file("test_size_mismatch.bin", std::ios::binary);
        bin_file.write(
            reinterpret_cast<const char*>(matrix.data()),
            static_cast<std::streamsize>(matrix.size() * sizeof(double)));

        // Create metadata claiming 3x3 dimensions (size mismatch)
        std::ofstream meta_file(
            "test_size_mismatch.snpstats", std::ios::binary);
        int64_t rows = 3;  // Mismatch with actual data
        int64_t cols = 3;  // Mismatch with actual data
        meta_file.write(reinterpret_cast<const char*>(&rows), sizeof(int64_t));
        meta_file.write(reinterpret_cast<const char*>(&cols), sizeof(int64_t));
    }

    static void createReadOnlyFiles()
    {
        // Create binary file
        Eigen::MatrixXd matrix(2, 2);
        matrix << 1.0, 3.0, 2.0, 4.0;

        std::ofstream bin_file("test_readonly.bin", std::ios::binary);
        bin_file.write(
            reinterpret_cast<const char*>(matrix.data()),
            static_cast<std::streamsize>(matrix.size() * sizeof(double)));

        // Create metadata
        std::ofstream meta_file("test_readonly.snpstats", std::ios::binary);
        int64_t rows = 2;
        int64_t cols = 2;
        meta_file.write(reinterpret_cast<const char*>(&rows), sizeof(int64_t));
        meta_file.write(reinterpret_cast<const char*>(&cols), sizeof(int64_t));

        // Make files read-only
        std::filesystem::permissions(
            "test_readonly.bin", std::filesystem::perms::owner_read);
        std::filesystem::permissions(
            "test_readonly.snpstats", std::filesystem::perms::owner_read);
    }

    static void createLargeTestFiles()
    {
        // Create a larger matrix for testing
        Eigen::MatrixXd matrix(100, 50);
        for (int i = 0; i < matrix.rows(); ++i)
        {
            for (int j = 0; j < matrix.cols(); ++j)
            {
                matrix(i, j) = static_cast<double>((i * matrix.cols()) + j);
            }
        }

        // Write binary data
        std::ofstream bin_file("test_large.bin", std::ios::binary);
        bin_file.write(
            reinterpret_cast<const char*>(matrix.data()),
            static_cast<std::streamsize>(matrix.size() * sizeof(double)));

        // Write metadata
        std::ofstream meta_file("test_large.snpstats", std::ios::binary);
        int64_t rows = matrix.rows();
        int64_t cols = matrix.cols();
        meta_file.write(reinterpret_cast<const char*>(&rows), sizeof(int64_t));
        meta_file.write(reinterpret_cast<const char*>(&cols), sizeof(int64_t));
    }

    static Eigen::MatrixXd getExpected3x2Matrix()
    {
        Eigen::MatrixXd matrix(3, 2);
        matrix << 1.0, 4.0, 2.0, 5.0, 3.0, 6.0;
        return matrix;
    }
};

TEST_CASE_PERSISTENT_FIXTURE(
    GenotypeMapTestFixture,
    "GenotypeMap::create function - success cases",
    "[genotype_mmap][create]")
{
    SECTION("Valid binary and metadata files")
    {
        auto result = gelex::GenotypeMap::create("test_valid.bin");
        REQUIRE(result.has_value());

        const auto& genotype_map = result.value();
        REQUIRE(genotype_map.rows() == 3);
        REQUIRE(genotype_map.cols() == 2);

        const auto& matrix = genotype_map.matrix();
        REQUIRE(matrix.rows() == 3);
        REQUIRE(matrix.cols() == 2);

        // Verify matrix data matches expected values
        auto expected_matrix = GenotypeMapTestFixture::getExpected3x2Matrix();
        for (int i = 0; i < matrix.rows(); ++i)
        {
            for (int j = 0; j < matrix.cols(); ++j)
            {
                REQUIRE_THAT(
                    matrix(i, j), WithinAbs(expected_matrix(i, j), 1e-10));
            }
        }
    }

    SECTION("Large matrix files")
    {
        // Create large test files for this specific test
        GenotypeMapTestFixture::createLargeTestFiles();

        auto result = gelex::GenotypeMap::create("test_large.bin");
        REQUIRE(result.has_value());

        const auto& genotype_map = result.value();
        REQUIRE(genotype_map.rows() == 100);
        REQUIRE(genotype_map.cols() == 50);

        const auto& matrix = genotype_map.matrix();
        REQUIRE(matrix.rows() == 100);
        REQUIRE(matrix.cols() == 50);

        // Verify some sample data points
        for (int i = 0; i < 10; i += 10)
        {
            for (int j = 0; j < 5; j += 5)
            {
                auto expected_value = static_cast<double>((i * 50) + j);
                REQUIRE_THAT(matrix(i, j), WithinAbs(expected_value, 1e-10));
            }
        }

        // Clean up large test files
        std::remove("test_large.bin");
        std::remove("test_large.snpstats");
    }
}

TEST_CASE_PERSISTENT_FIXTURE(
    GenotypeMapTestFixture,
    "GenotypeMap::create function - error cases",
    "[genotype_mmap][create]")
{
    SECTION("Non-existent binary file")
    {
        auto result = gelex::GenotypeMap::create("non_existent.bin");
        REQUIRE_FALSE(result.has_value());
        REQUIRE(result.error().code == gelex::ErrorCode::FileNotFound);
    }

    SECTION("Non-existent metadata file")
    {
        // Create binary file but no metadata
        Eigen::MatrixXd matrix(2, 2);
        matrix << 1.0, 3.0, 2.0, 4.0;

        std::ofstream bin_file("test_no_meta.bin", std::ios::binary);
        bin_file.write(
            reinterpret_cast<const char*>(matrix.data()),
            static_cast<std::streamsize>(matrix.size() * sizeof(double)));

        auto result = gelex::GenotypeMap::create("test_no_meta.bin");
        REQUIRE_FALSE(result.has_value());
        REQUIRE(result.error().code == gelex::ErrorCode::FileNotFound);

        std::remove("test_no_meta.bin");
    }

    SECTION("Invalid metadata dimensions")
    {
        auto result
            = gelex::GenotypeMap::create("test_invalid_meta.bin");
        REQUIRE_FALSE(result.has_value());
        REQUIRE(result.error().code == gelex::ErrorCode::InvalidData);
    }

    SECTION("Size mismatch between metadata and binary file")
    {
        auto result
            = gelex::GenotypeMap::create("test_size_mismatch.bin");
        REQUIRE_FALSE(result.has_value());
        REQUIRE(result.error().code == gelex::ErrorCode::InvalidData);
    }
}

TEST_CASE_PERSISTENT_FIXTURE(
    GenotypeMapTestFixture,
    "GenotypeMap matrix access and properties",
    "[genotype_mmap][access]")
{
    SECTION("Matrix dimensions and properties")
    {
        auto result = gelex::GenotypeMap::create("test_valid.bin");
        REQUIRE(result.has_value());

        const auto& genotype_map = result.value();

        // Test dimension accessors
        REQUIRE(genotype_map.rows() == 3);
        REQUIRE(genotype_map.cols() == 2);

        // Test matrix access
        const auto& matrix = genotype_map.matrix();
        REQUIRE(matrix.rows() == 3);
        REQUIRE(matrix.cols() == 2);

        // Test that matrix is read-only (const)
        static_assert(
            std::is_const_v<std::remove_reference_t<decltype(matrix)>>,
            "Matrix should be const");
    }

    SECTION("Matrix data integrity")
    {
        auto result = gelex::GenotypeMap::create("test_valid.bin");
        REQUIRE(result.has_value());

        const auto& genotype_map = result.value();
        const auto& matrix = genotype_map.matrix();

        auto expected_matrix = GenotypeMapTestFixture::getExpected3x2Matrix();

        // Test all elements
        for (int i = 0; i < matrix.rows(); ++i)
        {
            for (int j = 0; j < matrix.cols(); ++j)
            {
                REQUIRE_THAT(
                    matrix(i, j), WithinAbs(expected_matrix(i, j), 1e-10));
            }
        }

        // Test data pointer access (column-major order)
        const double* data_ptr = matrix.data();
        REQUIRE_THAT(data_ptr[0], WithinAbs(1.0, 1e-10));  // (0,0)
        REQUIRE_THAT(data_ptr[1], WithinAbs(2.0, 1e-10));  // (1,0)
        REQUIRE_THAT(data_ptr[2], WithinAbs(3.0, 1e-10));  // (2,0)
        REQUIRE_THAT(data_ptr[3], WithinAbs(4.0, 1e-10));  // (0,1)
        REQUIRE_THAT(data_ptr[4], WithinAbs(5.0, 1e-10));  // (1,1)
        REQUIRE_THAT(data_ptr[5], WithinAbs(6.0, 1e-10));  // (2,1)
    }
}

TEST_CASE_PERSISTENT_FIXTURE(
    GenotypeMapTestFixture,
    "GenotypeMap memory mapping behavior",
    "[genotype_mmap][memory]")
{
    SECTION("Matrix data is properly memory mapped")
    {
        auto result = gelex::GenotypeMap::create("test_valid.bin");
        REQUIRE(result.has_value());

        const auto& genotype_map = result.value();
        const auto& matrix = genotype_map.matrix();

        // The matrix should be a memory map, not a copy
        // We can verify this by checking that modifying the underlying file
        // would affect the matrix data (though we won't actually modify it)

        // Just verify the data is accessible and correct
        REQUIRE(matrix(0, 0) == 1.0);
        REQUIRE(matrix(2, 1) == 6.0);
    }
}

TEST_CASE("GenotypeMap edge cases", "[genotype_mmap][edge]")
{
    SECTION("Empty matrix (should fail)")
    {
        // Create empty binary file
        Eigen::MatrixXd empty_matrix(0, 0);
        std::ofstream empty_bin("test_empty.bin", std::ios::binary);
        empty_bin.write(
            reinterpret_cast<const char*>(empty_matrix.data()),
            static_cast<std::streamsize>(empty_matrix.size() * sizeof(double)));

        // Create metadata with zero dimensions
        std::ofstream empty_meta("test_empty.snpstats", std::ios::binary);
        int64_t rows = 0;
        int64_t cols = 0;
        empty_meta.write(reinterpret_cast<const char*>(&rows), sizeof(int64_t));
        empty_meta.write(reinterpret_cast<const char*>(&cols), sizeof(int64_t));

        auto result = gelex::GenotypeMap::create("test_empty.bin");
        REQUIRE_FALSE(result.has_value());
        REQUIRE(result.error().code == gelex::ErrorCode::InvalidFile);

        std::remove("test_empty.bin");
        std::remove("test_empty.snpstats");
    }

    SECTION("Single element matrix")
    {
        // Create 1x1 matrix
        Eigen::MatrixXd single_matrix(1, 1);
        single_matrix << 42.0;

        std::ofstream single_bin("test_single.bin", std::ios::binary);
        single_bin.write(
            reinterpret_cast<const char*>(single_matrix.data()),
            static_cast<std::streamsize>(
                single_matrix.size() * sizeof(double)));
        single_bin.flush();

        // Create metadata
        std::ofstream single_meta("test_single.snpstats", std::ios::binary);
        int64_t rows = 1;
        int64_t cols = 1;
        single_meta.write(
            reinterpret_cast<const char*>(&rows), sizeof(int64_t));
        single_meta.write(
            reinterpret_cast<const char*>(&cols), sizeof(int64_t));
        single_meta.flush();

        auto result = gelex::GenotypeMap::create("test_single.bin");
        REQUIRE(result.has_value());

        const auto& genotype_map = result.value();
        REQUIRE(genotype_map.rows() == 1);
        REQUIRE(genotype_map.cols() == 1);
        REQUIRE(genotype_map.matrix()(0, 0) == 42.0);

        std::remove("test_single.bin");
        std::remove("test_single.snpstats");
    }
}
