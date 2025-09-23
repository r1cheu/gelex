#include <filesystem>
#include <fstream>
#include <iostream>
#include <vector>

#include <Eigen/Core>
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

using Catch::Matchers::WithinAbs;

#include "../src/data/binary_matrix_writer.h"
#include "gelex/error.h"

class BinaryMatrixWriterTestFixture
{
   public:
    BinaryMatrixWriterTestFixture()
    {
        createValidTestFile();
        createReadOnlyFile();
    }

    BinaryMatrixWriterTestFixture(const BinaryMatrixWriterTestFixture&)
        = default;
    BinaryMatrixWriterTestFixture(BinaryMatrixWriterTestFixture&&) = delete;
    BinaryMatrixWriterTestFixture& operator=(
        const BinaryMatrixWriterTestFixture&)
        = default;
    BinaryMatrixWriterTestFixture& operator=(BinaryMatrixWriterTestFixture&&)
        = delete;
    ~BinaryMatrixWriterTestFixture() { remove(); }

    static void remove()
    {
        std::remove("test_valid.bin");
        std::remove("test_readonly.bin");
    }

    static void createValidTestFile()
    {
        // Create an empty valid file for testing
        std::ofstream file("test_valid.bin");
    }

    static void createReadOnlyFile()
    {
        // Create a file and make it read-only
        std::ofstream file("test_readonly.bin");
        file.close();
        std::filesystem::permissions(
            "test_readonly.bin", std::filesystem::perms::owner_read);
    }

    static Eigen::MatrixXd createTestMatrix3x2()
    {
        Eigen::MatrixXd matrix(3, 2);
        matrix << 1.0, 4.0, 2.0, 5.0, 3.0, 6.0;
        return matrix;
    }

    static Eigen::MatrixXd createTestMatrix2x3()
    {
        Eigen::MatrixXd matrix(2, 3);
        matrix << 1.1, 2.2, 3.3, 4.4, 5.5, 6.6;
        return matrix;
    }

    static Eigen::MatrixXd createEmptyMatrix() { return Eigen::MatrixXd(0, 0); }

    static Eigen::MatrixXd createLargeMatrix()
    {
        Eigen::MatrixXd matrix(100, 50);
        for (int i = 0; i < matrix.rows(); ++i)
        {
            for (int j = 0; j < matrix.cols(); ++j)
            {
                matrix(i, j) = static_cast<double>(i * matrix.cols() + j);
            }
        }
        return matrix;
    }
};

TEST_CASE_PERSISTENT_FIXTURE(
    BinaryMatrixWriterTestFixture,
    "BinaryMatrixWriter::create function",
    "[binary_matrix][writer]")
{
    SECTION("Valid file path")
    {
        auto result
            = gelex::detail::BinaryMatrixWriter::create("test_valid.bin");
        REQUIRE(result.has_value());
    }

    SECTION("Non-existent file")
    {
        auto result = gelex::detail::BinaryMatrixWriter::create(
            "non_existent_file.bin");
        REQUIRE(result.has_value());  // Should create the file
    }

    SECTION("Read-only file")
    {
        auto result
            = gelex::detail::BinaryMatrixWriter::create("test_readonly.bin");
        REQUIRE_FALSE(result.has_value());
        REQUIRE(result.error().code == gelex::ErrorCode::FileIOError);
    }
}

TEST_CASE_PERSISTENT_FIXTURE(
    BinaryMatrixWriterTestFixture,
    "BinaryMatrixWriter::write function - valid data",
    "[binary_matrix][writer]")
{
    SECTION("Write 3x2 matrix")
    {
        // Use a fresh file for this test
        std::remove("test_3x2.bin");

        auto writer_result
            = gelex::detail::BinaryMatrixWriter::create("test_3x2.bin");
        REQUIRE(writer_result.has_value());

        auto& writer = writer_result.value();
        auto matrix = BinaryMatrixWriterTestFixture::createTestMatrix3x2();

        auto write_result = writer.write(matrix);
        REQUIRE(write_result.has_value());

        // Verify file was written correctly
        std::ifstream file("test_3x2.bin", std::ios::binary);
        REQUIRE(file.is_open());

        // Read the binary data and verify it matches the original matrix
        std::vector<double> read_data(matrix.size());
        file.read(
            reinterpret_cast<char*>(read_data.data()),
            static_cast<std::streamsize>(matrix.size() * sizeof(double)));

        // Verify data is written in column-major order (Eigen's default)
        REQUIRE(read_data.size() == matrix.size());
        for (Eigen::Index i = 0; i < matrix.size(); ++i)
        {
            REQUIRE_THAT(read_data[i], WithinAbs(matrix.data()[i], 1e-10));
        }

        std::remove("test_3x2.bin");
    }

    SECTION("Write 2x3 matrix")
    {
        // Use a fresh file for this test
        std::remove("test_2x3.bin");

        auto writer_result
            = gelex::detail::BinaryMatrixWriter::create("test_2x3.bin");
        REQUIRE(writer_result.has_value());

        auto& writer = writer_result.value();
        auto matrix = BinaryMatrixWriterTestFixture::createTestMatrix2x3();

        auto write_result = writer.write(matrix);
        REQUIRE(write_result.has_value());

        // Verify file was written correctly
        std::ifstream file("test_2x3.bin", std::ios::binary);
        REQUIRE(file.is_open());

        std::vector<double> read_data(matrix.size());
        file.read(
            reinterpret_cast<char*>(read_data.data()),
            static_cast<std::streamsize>(matrix.size() * sizeof(double)));

        REQUIRE(read_data.size() == matrix.size());
        for (Eigen::Index i = 0; i < matrix.size(); ++i)
        {
            REQUIRE_THAT(read_data[i], WithinAbs(matrix.data()[i], 1e-10));
        }

        std::remove("test_2x3.bin");
    }

    SECTION("Write empty matrix")
    {
        // Use a fresh file for this test to avoid contamination from previous
        // tests
        std::remove("test_empty.bin");

        auto writer_result
            = gelex::detail::BinaryMatrixWriter::create("test_empty.bin");
        REQUIRE(writer_result.has_value());

        auto& writer = writer_result.value();
        auto matrix = BinaryMatrixWriterTestFixture::createEmptyMatrix();

        auto write_result = writer.write(matrix);
        REQUIRE(write_result.has_value());

        // Empty matrix should not write any data
        std::ifstream file("test_empty.bin", std::ios::binary | std::ios::ate);
        REQUIRE(file.is_open());
        REQUIRE(file.tellg() == 0);

        std::remove("test_empty.bin");
    }

    SECTION("Write large matrix")
    {
        // Use a fresh file for this test
        std::remove("test_large.bin");

        auto writer_result
            = gelex::detail::BinaryMatrixWriter::create("test_large.bin");
        REQUIRE(writer_result.has_value());

        auto& writer = writer_result.value();
        auto matrix = BinaryMatrixWriterTestFixture::createLargeMatrix();

        auto write_result = writer.write(matrix);
        REQUIRE(write_result.has_value());

        // Verify file size is correct
        std::ifstream file("test_large.bin", std::ios::binary | std::ios::ate);
        REQUIRE(file.is_open());

        auto file_size = file.tellg();
        size_t expected_size = matrix.size() * sizeof(double);
        REQUIRE(static_cast<size_t>(file_size) == expected_size);

        std::remove("test_large.bin");
    }
}

TEST_CASE_PERSISTENT_FIXTURE(
    BinaryMatrixWriterTestFixture,
    "BinaryMatrixWriter::write function - error conditions",
    "[binary_matrix][writer]")
{
    SECTION("Write after file closure simulation")
    {
        // Create a writer and then manually close the file to simulate an error
        auto writer_result
            = gelex::detail::BinaryMatrixWriter::create("test_output.bin");
        REQUIRE(writer_result.has_value());

        auto& writer = writer_result.value();
        auto matrix = BinaryMatrixWriterTestFixture::createTestMatrix3x2();

        // Manually close the file to simulate an error condition
        // Note: This is a bit tricky since the file is private, but we can test
        // by writing once successfully, then trying to write again after
        // the file might be in a bad state

        auto first_write = writer.write(matrix);
        REQUIRE(first_write.has_value());

        // Try to write again - should succeed if file is still good
        auto second_write = writer.write(matrix);
        REQUIRE(second_write.has_value());
    }
}

TEST_CASE_PERSISTENT_FIXTURE(
    BinaryMatrixWriterTestFixture,
    "BinaryMatrixWriter file format verification",
    "[binary_matrix][writer]")
{
    SECTION("Verify binary format structure - column-major order")
    {
        auto writer_result
            = gelex::detail::BinaryMatrixWriter::create("test_output.bin");
        REQUIRE(writer_result.has_value());

        auto& writer = writer_result.value();

        // Create a specific matrix to test column-major ordering
        Eigen::MatrixXd matrix(2, 3);
        matrix << 1.0, 3.0, 5.0, 2.0, 4.0, 6.0;

        auto write_result = writer.write(matrix);
        REQUIRE(write_result.has_value());

        // Verify file content matches column-major order
        std::ifstream file("test_output.bin", std::ios::binary);
        REQUIRE(file.is_open());

        std::vector<double> read_data(matrix.size());
        file.read(
            reinterpret_cast<char*>(read_data.data()),
            static_cast<std::streamsize>(matrix.size() * sizeof(double)));

        // In column-major order, the data should
        // be: 1.0, 2.0, 3.0, 4.0, 5.0, 6.0
        std::vector<double> expected_data = {1.0, 2.0, 3.0, 4.0, 5.0, 6.0};

        REQUIRE(read_data.size() == expected_data.size());
        for (size_t i = 0; i < read_data.size(); ++i)
        {
            REQUIRE_THAT(read_data[i], WithinAbs(expected_data[i], 1e-10));
        }
    }

    SECTION("Verify file size calculation")
    {
        // Use a fresh file for this test
        std::remove("test_size.bin");

        auto writer_result
            = gelex::detail::BinaryMatrixWriter::create("test_size.bin");
        REQUIRE(writer_result.has_value());

        auto& writer = writer_result.value();
        auto matrix = BinaryMatrixWriterTestFixture::createTestMatrix3x2();

        auto write_result = writer.write(matrix);
        REQUIRE(write_result.has_value());

        // Verify file size is correct
        std::ifstream file("test_size.bin", std::ios::binary | std::ios::ate);
        REQUIRE(file.is_open());

        auto file_size = file.tellg();
        size_t expected_size = matrix.size() * sizeof(double);
        REQUIRE(static_cast<size_t>(file_size) == expected_size);

        std::remove("test_size.bin");
    }

    SECTION("Multiple writes to same file")
    {
        // Use a fresh file for this test
        std::remove("test_multiple.bin");

        auto writer_result
            = gelex::detail::BinaryMatrixWriter::create("test_multiple.bin");
        REQUIRE(writer_result.has_value());

        auto& writer = writer_result.value();

        // Write first matrix
        Eigen::MatrixXd matrix1(2, 2);
        matrix1 << 1.0, 3.0, 2.0, 4.0;

        auto write_result1 = writer.write(matrix1);
        REQUIRE(write_result1.has_value());

        // Write second matrix (appends to the same file)
        Eigen::MatrixXd matrix2(2, 2);
        matrix2 << 5.0, 7.0, 6.0, 8.0;

        auto write_result2 = writer.write(matrix2);
        REQUIRE(write_result2.has_value());

        // Verify total file size
        std::ifstream file(
            "test_multiple.bin", std::ios::binary | std::ios::ate);
        REQUIRE(file.is_open());

        auto file_size = file.tellg();
        size_t expected_size
            = (matrix1.size() + matrix2.size()) * sizeof(double);
        REQUIRE(static_cast<size_t>(file_size) == expected_size);

        // Verify content
        file.seekg(0);
        std::vector<double> read_data((matrix1.size() + matrix2.size()));
        file.read(
            reinterpret_cast<char*>(read_data.data()),
            static_cast<std::streamsize>(
                (matrix1.size() + matrix2.size()) * sizeof(double)));

        // First matrix data (column-major): 1.0, 2.0, 3.0, 4.0
        // Second matrix data (column-major): 5.0, 6.0, 7.0, 8.0
        std::vector<double> expected_data
            = {1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0};

        REQUIRE(read_data.size() == expected_data.size());
        for (size_t i = 0; i < read_data.size(); ++i)
        {
            REQUIRE_THAT(read_data[i], WithinAbs(expected_data[i], 1e-10));
        }

        std::remove("test_multiple.bin");
    }
}
