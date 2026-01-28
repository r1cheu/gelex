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

#include <math.h>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <filesystem>
#include <fstream>
#include <limits>
#include <vector>

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_exception.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>

#include "file_fixture.h"
#include "gelex/exception.h"

#include "../src/data/binary_matrix_writer.h"

namespace fs = std::filesystem;

using namespace gelex::detail;  // NOLINT
using Catch::Matchers::EndsWith;
using Catch::Matchers::MessageMatches;
using gelex::test::FileFixture;

TEST_CASE(
    "BinaryMatrixWriter - Constructor and path access",
    "[data][binary_matrix]")
{
    FileFixture files;

    SECTION("Happy path - create writer with valid path")
    {
        auto file_path = files.generate_random_file_path(".bin");

        REQUIRE_NOTHROW(
            [&]()
            {
                BinaryMatrixWriter writer(file_path);
                REQUIRE(writer.path() == file_path);
            }());
    }
}

TEST_CASE("BinaryMatrixWriter - Write empty matrix", "[data][binary_matrix]")
{
    FileFixture files;

    SECTION("Happy path - write empty matrix (should skip)")
    {
        auto file_path = files.generate_random_file_path(".bin");

        Eigen::MatrixXd empty_matrix(0, 0);

        REQUIRE_NOTHROW(
            [&]()
            {
                BinaryMatrixWriter writer(file_path);
                writer.write(empty_matrix);
            }());

        REQUIRE(fs::exists(file_path));
        REQUIRE(fs::file_size(file_path) == 0);
    }

    SECTION("Happy path - write zero-sized matrix (0x0)")
    {
        auto file_path = files.generate_random_file_path(".bin");

        Eigen::MatrixXd zero_matrix(0, 0);

        REQUIRE_NOTHROW(
            [&]()
            {
                BinaryMatrixWriter writer(file_path);
                writer.write(zero_matrix);
            }());

        REQUIRE(fs::exists(file_path));
        REQUIRE(fs::file_size(file_path) == 0);
    }

    SECTION("Happy path - write matrix with zero rows (0xn)")
    {
        auto file_path = files.generate_random_file_path(".bin");

        Eigen::MatrixXd zero_rows_matrix(0, 5);

        REQUIRE_NOTHROW(
            [&]()
            {
                BinaryMatrixWriter writer(file_path);
                writer.write(zero_rows_matrix);
            }());

        REQUIRE(fs::exists(file_path));
        REQUIRE(fs::file_size(file_path) == 0);
    }

    SECTION("Happy path - write matrix with zero columns (nx0)")
    {
        auto file_path = files.generate_random_file_path(".bin");

        Eigen::MatrixXd zero_cols_matrix(5, 0);

        REQUIRE_NOTHROW(
            [&]()
            {
                BinaryMatrixWriter writer(file_path);
                writer.write(zero_cols_matrix);
            }());

        REQUIRE(fs::exists(file_path));
        REQUIRE(fs::file_size(file_path) == 0);
    }
}

TEST_CASE(
    "BinaryMatrixWriter - Write non-empty matrix",
    "[data][binary_matrix]")
{
    FileFixture files;

    SECTION(
        "Happy path - write basic matrix with positive, negative, and zero "
        "values")
    {
        auto file_path = files.generate_random_file_path(".bin");

        Eigen::MatrixXd matrix(3, 3);
        matrix << 1.5, -2.0, 3.0, -4.5, 0.0, 6.0, 7.0, -8.5, 9.0;

        REQUIRE_NOTHROW(
            [&]()
            {
                BinaryMatrixWriter writer(file_path);
                writer.write(matrix);
            }());

        REQUIRE(fs::exists(file_path));
        REQUIRE(fs::file_size(file_path) == matrix.size() * sizeof(double));

        std::ifstream file(file_path, std::ios::binary);
        REQUIRE(file.is_open());

        std::vector<double> read_data(matrix.size());
        file.read(
            reinterpret_cast<char*>(read_data.data()),
            static_cast<std::streamsize>(matrix.size() * sizeof(double)));

        REQUIRE(read_data.size() == matrix.size());
        for (int i = 0; i < matrix.size(); ++i)
        {
            REQUIRE(read_data[i] == matrix.data()[i]);
        }
    }

    SECTION("Happy path - write 1x1 matrix (scalar)")
    {
        auto file_path = files.generate_random_file_path(".bin");

        Eigen::MatrixXd matrix(1, 1);
        matrix << 42.0;

        REQUIRE_NOTHROW(
            [&]()
            {
                BinaryMatrixWriter writer(file_path);
                writer.write(matrix);
            }());

        REQUIRE(fs::exists(file_path));
        REQUIRE(fs::file_size(file_path) == sizeof(double));

        std::ifstream file(file_path, std::ios::binary);
        REQUIRE(file.is_open());

        double read_value = {};
        file.read(reinterpret_cast<char*>(&read_value), sizeof(double));
        REQUIRE(read_value == 42.0);
    }

    SECTION("Happy path - write medium matrix 10x10")
    {
        auto file_path = files.generate_random_file_path(".bin");

        const int rows = 10;
        const int cols = 10;
        Eigen::MatrixXd matrix = Eigen::MatrixXd::Random(rows, cols);

        REQUIRE_NOTHROW(
            [&]()
            {
                BinaryMatrixWriter writer(file_path);
                writer.write(matrix);
            }());

        REQUIRE(fs::exists(file_path));
        REQUIRE(fs::file_size(file_path) == matrix.size() * sizeof(double));

        std::ifstream file(file_path, std::ios::binary);
        REQUIRE(file.is_open());

        std::vector<double> read_data(matrix.size());
        file.read(
            reinterpret_cast<char*>(read_data.data()),
            static_cast<std::streamsize>(matrix.size() * sizeof(double)));

        REQUIRE(read_data.size() == matrix.size());
        for (int i = 0; i < matrix.size(); ++i)
        {
            REQUIRE(read_data[i] == matrix.data()[i]);
        }
    }

    SECTION("Happy path - write matrix with special values (inf, nan)")
    {
        auto file_path = files.generate_random_file_path(".bin");

        Eigen::MatrixXd matrix(3, 3);
        const double inf = std::numeric_limits<double>::infinity();
        const double neg_inf = -std::numeric_limits<double>::infinity();
        const double nan = std::numeric_limits<double>::quiet_NaN();

        matrix << 1.0, 2.0, 3.0, inf, neg_inf, 4.0, 5.0, nan, 6.0;

        REQUIRE_NOTHROW(
            [&]()
            {
                BinaryMatrixWriter writer(file_path);
                writer.write(matrix);
            }());

        REQUIRE(fs::exists(file_path));
        REQUIRE(fs::file_size(file_path) == matrix.size() * sizeof(double));

        std::ifstream file(file_path, std::ios::binary);
        REQUIRE(file.is_open());

        std::vector<double> read_data(matrix.size());
        file.read(
            reinterpret_cast<char*>(read_data.data()),
            static_cast<std::streamsize>(matrix.size() * sizeof(double)));

        REQUIRE(read_data.size() == matrix.size());

        REQUIRE(read_data[0] == 1.0);
        REQUIRE((std::isinf(read_data[1]) && read_data[1] > 0));
        REQUIRE(read_data[2] == 5.0);

        // 列1: 2.0, -inf, nan
        REQUIRE(read_data[3] == 2.0);
        REQUIRE((std::isinf(read_data[4]) && read_data[4] < 0));
        REQUIRE(std::isnan(read_data[5]));

        // 列2: 3.0, 4.0, 6.0
        REQUIRE(read_data[6] == 3.0);
        REQUIRE(read_data[7] == 4.0);
        REQUIRE(read_data[8] == 6.0);
    }

    SECTION("Happy path - multiple writes to same writer")
    {
        auto file_path = files.generate_random_file_path(".bin");

        Eigen::MatrixXd matrix1(2, 2);
        matrix1 << 1.0, 2.0, 3.0, 4.0;

        Eigen::MatrixXd matrix2(2, 2);
        matrix2 << 5.0, 6.0, 7.0, 8.0;

        REQUIRE_NOTHROW(
            [&]()
            {
                BinaryMatrixWriter writer(file_path);
                writer.write(matrix1);
                writer.write(matrix2);
            }());

        REQUIRE(fs::exists(file_path));
        REQUIRE(
            fs::file_size(file_path)
            == (matrix1.size() + matrix2.size()) * sizeof(double));

        std::ifstream file(file_path, std::ios::binary);
        REQUIRE(file.is_open());

        std::vector<double> read_data(matrix1.size() + matrix2.size());
        file.read(
            reinterpret_cast<char*>(read_data.data()),
            static_cast<std::streamsize>(read_data.size() * sizeof(double)));

        for (int i = 0; i < matrix1.size(); ++i)
        {
            REQUIRE(read_data[i] == matrix1.data()[i]);
        }

        for (int i = 0; i < matrix2.size(); ++i)
        {
            REQUIRE(read_data[matrix1.size() + i] == matrix2.data()[i]);
        }
    }
}

TEST_CASE(
    "BinaryMatrixWriter - Buffer size verification",
    "[data][binary_matrix]")
{
    FileFixture files;

    SECTION("Happy path - verify default buffer size constant")
    {
        REQUIRE(
            BinaryMatrixWriter::kDefaultBufferSize
            == static_cast<size_t>(64 * 1024));
    }

    SECTION("Happy path - write with default buffer should work")
    {
        auto file_path = files.generate_random_file_path(".bin");

        Eigen::MatrixXd matrix(10, 10);
        matrix = Eigen::MatrixXd::Random(10, 10);

        REQUIRE_NOTHROW(
            [&]()
            {
                BinaryMatrixWriter writer(file_path);
                writer.write(matrix);
            }());

        REQUIRE(fs::exists(file_path));
        REQUIRE(fs::file_size(file_path) == matrix.size() * sizeof(double));
    }
}
