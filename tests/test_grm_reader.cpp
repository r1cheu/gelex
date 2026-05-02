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

#include <cmath>
#include <cstddef>
#include <filesystem>
#include <string>
#include <vector>

#include <Eigen/Core>
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_exception.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>

#include "file_fixture.h"
#include "gelex/data/dataframe/index.h"
#include "gelex/data/grm/detail/grm_reader.h"
#include "gelex/data/grm/grm_writer.h"
#include "gelex/exception.h"
#include "gelex/types/sample_id.h"

namespace fs = std::filesystem;
namespace df = gelex::dataframe;

using namespace gelex::grm::detail;  // NOLINT
using Catch::Matchers::ContainsSubstring;
using Catch::Matchers::MessageMatches;
using gelex::test::FileFixture;

// Helper class to create GRM test files (bin + id)
class GrmFileFixture
{
   public:
    explicit GrmFileFixture(FileFixture& files)
        : files_(files), prefix_(files.generate_random_file_path(""))
    {
    }

    // Create GRM files from matrix and IDs
    auto create(
        const Eigen::MatrixXd& matrix,
        const std::vector<std::string>& ids) -> void
    {
        // Write binary file
        auto bin_path = fs::path(prefix_.string() + ".bin");
        {
            gelex::GrmBinWriter writer(bin_path);
            writer.write(matrix);
        }

        gelex::write_grm_ids(prefix_.string() + ".id", ids);
    }

    // Create only ID file (for testing missing bin file)
    auto create_id_only(const std::vector<std::string>& ids) -> void
    {
        gelex::write_grm_ids(prefix_.string() + ".id", ids);
    }

    // Create only bin file (for testing missing id file)
    auto create_bin_only(const Eigen::MatrixXd& matrix) -> void
    {
        auto bin_path = fs::path(prefix_.string() + ".bin");
        gelex::GrmBinWriter writer(bin_path);
        writer.write(matrix);
    }

    [[nodiscard]] auto prefix() const -> const fs::path& { return prefix_; }

   private:
    FileFixture& files_;
    fs::path prefix_;
};

// Helper to create a symmetric matrix
auto make_symmetric_matrix(Eigen::Index n) -> Eigen::MatrixXd
{
    Eigen::MatrixXd matrix = Eigen::MatrixXd::Random(n, n);
    return (matrix + matrix.transpose()) / 2.0;
}

// Helper to create canonical sample IDs
auto make_sample_ids(Eigen::Index n) -> std::vector<std::string>
{
    std::vector<std::string> ids;
    ids.reserve(static_cast<size_t>(n));
    for (Eigen::Index i = 0; i < n; ++i)
    {
        ids.push_back(
            gelex::make_sample_id(
                "FAM" + std::to_string(i), "IND" + std::to_string(i)));
    }
    return ids;
}

// ============================================================================
// Constructor tests
// ============================================================================

TEST_CASE(
    "GrmReader - Constructor with valid prefix",
    "[grm_reader][construction]")
{
    FileFixture files;
    GrmFileFixture grm_files(files);

    SECTION("Happy path - construct with valid files")
    {
        const Eigen::Index n = 3;
        auto matrix = make_symmetric_matrix(n);
        auto ids = make_sample_ids(n);

        grm_files.create(matrix, ids);

        REQUIRE_NOTHROW([&]() { GrmReader reader(grm_files.prefix()); }());
    }
}

TEST_CASE(
    "GrmReader - Constructor with missing files",
    "[grm_reader][construction][exception]")
{
    FileFixture files;
    GrmFileFixture grm_files(files);

    SECTION("Exception - missing .grm.bin file")
    {
        std::vector<std::string> ids
            = {gelex::make_sample_id("FAM1", "IND1"),
               gelex::make_sample_id("FAM2", "IND2")};
        grm_files.create_id_only(ids);

        REQUIRE_THROWS_AS(GrmReader(grm_files.prefix()), gelex::GelexException);
    }

    SECTION("Exception - missing .grm.id file")
    {
        Eigen::MatrixXd matrix(2, 2);
        matrix << 1.0, 0.5, 0.5, 1.0;

        grm_files.create_bin_only(matrix);

        REQUIRE_THROWS_AS(GrmReader(grm_files.prefix()), gelex::GelexException);
    }
}

TEST_CASE(
    "GrmReader - Constructor with size mismatch",
    "[grm_reader][construction][exception]")
{
    FileFixture files;
    GrmFileFixture grm_files(files);

    SECTION("Exception - bin file size doesn't match id count")
    {
        // Create a 3x3 matrix but only 2 IDs
        Eigen::MatrixXd matrix(3, 3);
        matrix << 1.0, 0.5, 0.3, 0.5, 1.0, 0.4, 0.3, 0.4, 1.0;

        std::vector<std::string> ids
            = {gelex::make_sample_id("FAM1", "IND1"),
               gelex::make_sample_id("FAM2", "IND2")};

        grm_files.create(matrix, ids);

        REQUIRE_THROWS_AS(GrmReader(grm_files.prefix()), gelex::GelexException);
    }

    SECTION("Exception - error message contains expected and actual size")
    {
        Eigen::MatrixXd matrix(4, 4);
        matrix.setIdentity();

        std::vector<std::string> ids
            = {gelex::make_sample_id("FAM1", "IND1"),
               gelex::make_sample_id("FAM2", "IND2")};

        grm_files.create(matrix, ids);

        REQUIRE_THROWS_MATCHES(
            GrmReader(grm_files.prefix()),
            gelex::GelexException,
            MessageMatches(ContainsSubstring("size mismatch")));
    }
}

// ============================================================================
// Accessor tests
// ============================================================================

TEST_CASE("GrmReader - sample_ids accessor", "[grm_reader][accessor]")
{
    FileFixture files;
    GrmFileFixture grm_files(files);

    SECTION("Happy path - verify loaded sample IDs")
    {
        const Eigen::Index n = 3;
        auto matrix = make_symmetric_matrix(n);
        std::vector<std::string> ids
            = {gelex::make_sample_id("FAM1", "IND1"),
               gelex::make_sample_id("FAM2", "IND2"),
               gelex::make_sample_id("FAM3", "IND3")};

        grm_files.create(matrix, ids);
        GrmReader reader(grm_files.prefix());

        const auto& idx = reader.sample_index();
        REQUIRE(idx.size() == 3);
        REQUIRE(idx.keys()[0] == gelex::make_sample_id("FAM1", "IND1"));
        REQUIRE(idx.keys()[1] == gelex::make_sample_id("FAM2", "IND2"));
        REQUIRE(idx.keys()[2] == gelex::make_sample_id("FAM3", "IND3"));
    }

    SECTION("Happy path - IDs with underscores in FID and IID are preserved")
    {
        const Eigen::Index n = 2;
        auto matrix = make_symmetric_matrix(n);
        std::vector<std::string> ids
            = {gelex::make_sample_id("FAM_1", "IND_1"),
               gelex::make_sample_id("FAM_2", "IND_2")};

        grm_files.create(matrix, ids);
        GrmReader reader(grm_files.prefix());

        const auto& idx = reader.sample_index();
        REQUIRE(idx.size() == 2);
        REQUIRE(idx.keys()[0] == gelex::make_sample_id("FAM_1", "IND_1"));
        REQUIRE(idx.keys()[1] == gelex::make_sample_id("FAM_2", "IND_2"));
    }
}

TEST_CASE("GrmReader - num_samples accessor", "[grm_reader][accessor]")
{
    FileFixture files;
    GrmFileFixture grm_files(files);

    SECTION("Happy path - verify num_samples for small matrix")
    {
        const Eigen::Index n = 5;
        auto matrix = make_symmetric_matrix(n);
        auto ids = make_sample_ids(n);

        grm_files.create(matrix, ids);
        GrmReader reader(grm_files.prefix());

        REQUIRE(reader.num_samples() == n);
    }

    SECTION("Happy path - verify num_samples for larger matrix")
    {
        const Eigen::Index n = 50;
        auto matrix = make_symmetric_matrix(n);
        auto ids = make_sample_ids(n);

        grm_files.create(matrix, ids);
        GrmReader reader(grm_files.prefix());

        REQUIRE(reader.num_samples() == n);
    }
}

// ============================================================================
// load() - Complete matrix loading tests
// ============================================================================

TEST_CASE("GrmReader - Load complete 3x3 GRM", "[grm_reader][load]")
{
    FileFixture files;
    GrmFileFixture grm_files(files);

    SECTION("Happy path - load and verify 3x3 matrix values")
    {
        // Create known symmetric matrix
        // clang-format off
        Eigen::MatrixXd original(3, 3);
        original << 1.0, 0.5, 0.3,
                    0.5, 1.0, 0.4,
                    0.3, 0.4, 1.0;
        // clang-format on

        auto ids = make_sample_ids(3);
        grm_files.create(original, ids);

        GrmReader reader(grm_files.prefix());
        auto loaded = reader.load();

        REQUIRE(loaded.rows() == 3);
        REQUIRE(loaded.cols() == 3);

        // Verify all elements (accounting for float precision loss)
        for (Eigen::Index i = 0; i < 3; ++i)
        {
            for (Eigen::Index j = 0; j < 3; ++j)
            {
                REQUIRE(
                    loaded(i, j)
                    == static_cast<double>(static_cast<float>(original(i, j))));
            }
        }
    }
}

TEST_CASE("GrmReader - Load complete 10x10 GRM", "[grm_reader][load]")
{
    FileFixture files;
    GrmFileFixture grm_files(files);

    SECTION("Happy path - load and verify 10x10 matrix (unnormalized)")
    {
        const Eigen::Index n = 10;
        auto original = make_symmetric_matrix(n);
        auto ids = make_sample_ids(n);

        grm_files.create(original, ids);

        GrmReader reader(grm_files.prefix());
        auto loaded = reader.load_unnormalized();

        REQUIRE(loaded.rows() == n);
        REQUIRE(loaded.cols() == n);

        // Verify all elements in unnormalized matrix
        for (Eigen::Index i = 0; i < n; ++i)
        {
            for (Eigen::Index j = 0; j < n; ++j)
            {
                auto expected
                    = static_cast<double>(static_cast<float>(original(i, j)));
                REQUIRE(loaded(i, j) == expected);
            }
        }
    }
}

TEST_CASE("GrmReader - Verify loaded matrix is symmetric", "[grm_reader][load]")
{
    FileFixture files;
    GrmFileFixture grm_files(files);

    SECTION("Happy path - loaded matrix should be symmetric")
    {
        const Eigen::Index n = 5;
        auto original = make_symmetric_matrix(n);
        auto ids = make_sample_ids(n);

        grm_files.create(original, ids);

        GrmReader reader(grm_files.prefix());
        auto loaded = reader.load();

        // Verify symmetry: M(i,j) == M(j,i)
        for (Eigen::Index i = 0; i < n; ++i)
        {
            for (Eigen::Index j = 0; j < i; ++j)
            {
                REQUIRE(loaded(i, j) == loaded(j, i));
            }
        }
    }
}

TEST_CASE(
    "GrmReader - Verify numerical precision (double-float-double)",
    "[grm_reader][load][numerical]")
{
    FileFixture files;
    GrmFileFixture grm_files(files);

    SECTION("Happy path - verify precision loss from float conversion")
    {
        // Use values that lose precision when converted to float
        Eigen::MatrixXd original(2, 2);
        original << 1.23456789012345, 0.98765432109876, 0.98765432109876,
            0.00000012345678;

        auto ids = make_sample_ids(2);
        grm_files.create(original, ids);

        GrmReader reader(grm_files.prefix());
        auto loaded = reader.load_unnormalized();

        // Values should match float precision
        REQUIRE(
            loaded(0, 0)
            == static_cast<double>(static_cast<float>(original(0, 0))));
        REQUIRE(
            loaded(1, 1)
            == static_cast<double>(static_cast<float>(original(1, 1))));

        // Verify precision is indeed reduced (not equal to original double)
        REQUIRE(loaded(0, 0) != original(0, 0));
    }
}

// ============================================================================
// load(id_map) - Filtered/reordered loading tests
// ============================================================================

TEST_CASE("GrmReader - Load with subset of IDs", "[grm_reader][load_filtered]")
{
    FileFixture files;
    GrmFileFixture grm_files(files);

    SECTION("Happy path - load subset (2 of 4 samples)")
    {
        // Create 4x4 matrix with distinct values
        Eigen::MatrixXd original(4, 4);
        // clang-format off
        original << 1.0, 0.1, 0.2, 0.3,
                    0.1, 2.0, 0.4, 0.5,
                    0.2, 0.4, 3.0, 0.6,
                    0.3, 0.5, 0.6, 4.0;
        // clang-format on

        std::vector<std::string> ids
            = {gelex::make_sample_id("FAM0", "IND0"),
               gelex::make_sample_id("FAM1", "IND1"),
               gelex::make_sample_id("FAM2", "IND2"),
               gelex::make_sample_id("FAM3", "IND3")};

        grm_files.create(original, ids);

        GrmReader reader(grm_files.prefix());

        // Load only samples 1 and 3, mapping to indices 0 and 1
        df::Index<std::string> sample_index(
            std::vector<std::string>{
                gelex::make_sample_id("FAM1", "IND1"),
                gelex::make_sample_id("FAM3", "IND3")});

        auto loaded = reader.load_unnormalized(sample_index);

        REQUIRE(loaded.rows() == 2);
        REQUIRE(loaded.cols() == 2);

        // Verify values (using float conversion)
        auto to_float = [](double v)
        { return static_cast<double>(static_cast<float>(v)); };

        // (0,0) should be original(1,1) = 2.0
        REQUIRE(loaded(0, 0) == to_float(2.0));
        // (1,1) should be original(3,3) = 4.0
        REQUIRE(loaded(1, 1) == to_float(4.0));
        // (0,1) should be original(1,3) = 0.5
        REQUIRE(loaded(0, 1) == to_float(0.5));
        // (1,0) should be original(3,1) = 0.5
        REQUIRE(loaded(1, 0) == to_float(0.5));
    }
}

TEST_CASE("GrmReader - Load with reordered IDs", "[grm_reader][load_filtered]")
{
    FileFixture files;
    GrmFileFixture grm_files(files);

    SECTION("Happy path - reverse order of samples")
    {
        // Create 3x3 matrix
        Eigen::MatrixXd original(3, 3);
        // clang-format off
        original << 1.0, 0.1, 0.2,
                    0.1, 2.0, 0.3,
                    0.2, 0.3, 3.0;
        // clang-format on

        std::vector<std::string> ids
            = {gelex::make_sample_id("FAM0", "IND0"),
               gelex::make_sample_id("FAM1", "IND1"),
               gelex::make_sample_id("FAM2", "IND2")};

        grm_files.create(original, ids);

        GrmReader reader(grm_files.prefix());

        // Reverse the order: original[2]->0, original[1]->1, original[0]->2
        df::Index<std::string> sample_index(
            std::vector<std::string>{
                gelex::make_sample_id("FAM2", "IND2"),
                gelex::make_sample_id("FAM1", "IND1"),
                gelex::make_sample_id("FAM0", "IND0")});

        auto loaded = reader.load_unnormalized(sample_index);

        REQUIRE(loaded.rows() == 3);
        REQUIRE(loaded.cols() == 3);

        auto to_float = [](double v)
        { return static_cast<double>(static_cast<float>(v)); };

        // Diagonal should be reversed: 3.0, 2.0, 1.0
        REQUIRE(loaded(0, 0) == to_float(3.0));
        REQUIRE(loaded(1, 1) == to_float(2.0));
        REQUIRE(loaded(2, 2) == to_float(1.0));

        // Off-diagonal: loaded(0,2) = original(2,0) = 0.2
        REQUIRE(loaded(0, 2) == to_float(0.2));
        REQUIRE(loaded(2, 0) == to_float(0.2));
    }
}

TEST_CASE("GrmReader - Load with empty id_map", "[grm_reader][load_filtered]")
{
    FileFixture files;
    GrmFileFixture grm_files(files);

    SECTION("Happy path - empty id_map returns empty matrix")
    {
        const Eigen::Index n = 3;
        auto original = make_symmetric_matrix(n);
        auto ids = make_sample_ids(n);

        grm_files.create(original, ids);

        GrmReader reader(grm_files.prefix());

        df::Index<std::string> empty_index;
        auto loaded = reader.load(empty_index);

        REQUIRE(loaded.rows() == 0);
        REQUIRE(loaded.cols() == 0);
    }
}

TEST_CASE(
    "GrmReader - Load with invalid ID throws exception",
    "[grm_reader][load_filtered][exception]")
{
    FileFixture files;
    GrmFileFixture grm_files(files);

    SECTION("Exception - ID not found in file")
    {
        const Eigen::Index n = 3;
        auto original = make_symmetric_matrix(n);
        auto ids = make_sample_ids(n);

        grm_files.create(original, ids);

        GrmReader reader(grm_files.prefix());

        df::Index<std::string> sample_index(
            std::vector<std::string>{
                gelex::make_sample_id("FAM0", "IND0"),
                gelex::make_sample_id("NONEXISTENT", "ID")});

        REQUIRE_THROWS_AS(reader.load(sample_index), gelex::GelexException);
    }

    SECTION("Exception - error message contains the invalid ID")
    {
        const Eigen::Index n = 2;
        auto original = make_symmetric_matrix(n);
        auto ids = make_sample_ids(n);

        grm_files.create(original, ids);

        GrmReader reader(grm_files.prefix());

        df::Index<std::string> sample_index(
            std::vector<std::string>{
                gelex::make_sample_id("MISSING", "SAMPLE")});

        REQUIRE_THROWS_MATCHES(
            reader.load(sample_index),
            gelex::GelexException,
            MessageMatches(ContainsSubstring("MISSING")));
    }
}

// ============================================================================
// ID parsing tests
// ============================================================================

TEST_CASE("GrmReader - ID parsing from file", "[grm_reader][id_parsing]")
{
    FileFixture files;
    GrmFileFixture grm_files(files);

    SECTION("Happy path - canonical IDs round-trip through .id file")
    {
        const Eigen::Index n = 2;
        auto matrix = make_symmetric_matrix(n);
        std::vector<std::string> ids
            = {gelex::make_sample_id("SAMPLE", "1"),
               gelex::make_sample_id("SAMPLE", "2")};

        grm_files.create(matrix, ids);
        GrmReader reader(grm_files.prefix());

        const auto& idx = reader.sample_index();
        REQUIRE(idx.keys()[0] == gelex::make_sample_id("SAMPLE", "1"));
        REQUIRE(idx.keys()[1] == gelex::make_sample_id("SAMPLE", "2"));
    }
}

// ============================================================================
// Round-trip verification tests
// ============================================================================

TEST_CASE(
    "GrmReader - Round-trip write and load verification",
    "[grm_reader][roundtrip]")
{
    FileFixture files;
    GrmFileFixture grm_files(files);

    SECTION("Happy path - verify data integrity through write/load cycle")
    {
        const Eigen::Index n = 20;
        auto original = make_symmetric_matrix(n);
        auto ids = make_sample_ids(n);

        grm_files.create(original, ids);

        GrmReader reader(grm_files.prefix());
        auto loaded = reader.load_unnormalized();

        // Verify dimensions
        REQUIRE(loaded.rows() == n);
        REQUIRE(loaded.cols() == n);

        // Verify all values match (with float precision)
        for (Eigen::Index i = 0; i < n; ++i)
        {
            for (Eigen::Index j = 0; j < n; ++j)
            {
                auto expected
                    = static_cast<double>(static_cast<float>(original(i, j)));
                REQUIRE(loaded(i, j) == expected);
            }
        }
    }

    SECTION("Happy path - larger matrix round-trip")
    {
        const Eigen::Index n = 100;
        auto original = make_symmetric_matrix(n);
        auto ids = make_sample_ids(n);

        grm_files.create(original, ids);

        GrmReader reader(grm_files.prefix());
        auto loaded = reader.load_unnormalized();

        REQUIRE(loaded.rows() == n);
        REQUIRE(loaded.cols() == n);

        // Spot check diagonal and some off-diagonal elements
        for (Eigen::Index i = 0; i < n; i += 10)
        {
            auto expected
                = static_cast<double>(static_cast<float>(original(i, i)));
            REQUIRE(loaded(i, i) == expected);
        }
    }
}

// ============================================================================
// Special values tests
// ============================================================================

TEST_CASE(
    "GrmReader - Load matrix with special values",
    "[grm_reader][special]")
{
    FileFixture files;
    GrmFileFixture grm_files(files);

    SECTION("Happy path - matrix with inf values")
    {
        const auto inf = std::numeric_limits<double>::infinity();

        Eigen::MatrixXd original(2, 2);
        original << inf, 0.5, 0.5, 1.0;

        auto ids = make_sample_ids(2);
        grm_files.create(original, ids);

        GrmReader reader(grm_files.prefix());
        auto loaded = reader.load_unnormalized();

        REQUIRE(std::isinf(loaded(0, 0)));
        REQUIRE(loaded(0, 0) > 0);
    }

    SECTION("Happy path - matrix with NaN values")
    {
        const auto nan = std::numeric_limits<double>::quiet_NaN();

        Eigen::MatrixXd original(2, 2);
        original << nan, 0.5, 0.5, 1.0;

        auto ids = make_sample_ids(2);
        grm_files.create(original, ids);

        GrmReader reader(grm_files.prefix());
        auto loaded = reader.load_unnormalized();

        REQUIRE(std::isnan(loaded(0, 0)));
    }
}

// ============================================================================
// Normalization tests (trace-based)
// ============================================================================

TEST_CASE(
    "GrmReader - Load with trace-based normalization",
    "[grm_reader][normalize]")
{
    FileFixture files;
    GrmFileFixture grm_files(files);

    SECTION("Happy path - automatic normalization using trace/n")
    {
        // Create unnormalized matrix with known trace
        // trace = 2.0 + 4.0 + 6.0 = 12.0, n = 3
        // denominator = trace/n = 12.0/3 = 4.0
        Eigen::MatrixXd unnormalized(3, 3);
        unnormalized << 2.0, 1.0, 0.6, 1.0, 4.0, 0.8, 0.6, 0.8, 6.0;

        auto ids = make_sample_ids(3);
        grm_files.create(unnormalized, ids);

        GrmReader reader(grm_files.prefix());
        auto loaded = reader.load();

        // Verify automatic normalization (accounting for float precision)
        auto to_float = [](double v)
        { return static_cast<double>(static_cast<float>(v)); };

        // denominator = trace/n = (2.0 + 4.0 + 6.0) / 3 = 4.0
        REQUIRE(loaded(0, 0) == to_float(2.0 / 4.0));  // 2.0 / 4.0 = 0.5
        REQUIRE(loaded(1, 1) == to_float(4.0 / 4.0));  // 4.0 / 4.0 = 1.0
        REQUIRE(loaded(2, 2) == to_float(6.0 / 4.0));  // 6.0 / 4.0 = 1.5
    }
}
