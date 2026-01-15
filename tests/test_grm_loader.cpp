#include <cmath>
#include <cstddef>
#include <filesystem>
#include <string>
#include <unordered_map>
#include <vector>

#include <Eigen/Core>
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_exception.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>

#include "../src/data/grm_bin_writer.h"
#include "../src/data/grm_id_writer.h"
#include "../src/data/grm_loader.h"
#include "file_fixture.h"
#include "gelex/exception.h"

namespace fs = std::filesystem;

using namespace gelex::detail;  // NOLINT
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
        auto bin_path = fs::path(prefix_.string() + ".grm.bin");
        {
            GrmBinWriter writer(bin_path);
            writer.write(matrix);
        }

        // Write ID file
        auto id_path = fs::path(prefix_.string() + ".grm.id");
        {
            GrmIdWriter writer(id_path);
            writer.write(ids);
        }
    }

    // Create only ID file (for testing missing bin file)
    auto create_id_only(const std::vector<std::string>& ids) -> void
    {
        auto id_path = fs::path(prefix_.string() + ".grm.id");
        GrmIdWriter writer(id_path);
        writer.write(ids);
    }

    // Create only bin file (for testing missing id file)
    auto create_bin_only(const Eigen::MatrixXd& matrix) -> void
    {
        auto bin_path = fs::path(prefix_.string() + ".grm.bin");
        GrmBinWriter writer(bin_path);
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

// Helper to create IDs in FID_IID format
auto make_sample_ids(Eigen::Index n) -> std::vector<std::string>
{
    std::vector<std::string> ids;
    ids.reserve(static_cast<size_t>(n));
    for (Eigen::Index i = 0; i < n; ++i)
    {
        ids.push_back("FAM" + std::to_string(i) + "_IND" + std::to_string(i));
    }
    return ids;
}

// ============================================================================
// Constructor tests
// ============================================================================

TEST_CASE(
    "GrmLoader - Constructor with valid prefix",
    "[grm_loader][construction]")
{
    FileFixture files;
    GrmFileFixture grm_files(files);

    SECTION("Happy path - construct with valid files")
    {
        const Eigen::Index n = 3;
        auto matrix = make_symmetric_matrix(n);
        auto ids = make_sample_ids(n);

        grm_files.create(matrix, ids);

        REQUIRE_NOTHROW([&]() { GrmLoader loader(grm_files.prefix()); }());
    }
}

TEST_CASE(
    "GrmLoader - Constructor with missing files",
    "[grm_loader][construction][exception]")
{
    FileFixture files;
    GrmFileFixture grm_files(files);

    SECTION("Exception - missing .grm.bin file")
    {
        std::vector<std::string> ids = {"FAM1_IND1", "FAM2_IND2"};
        grm_files.create_id_only(ids);

        REQUIRE_THROWS_AS(
            GrmLoader(grm_files.prefix()), gelex::FileOpenException);
    }

    SECTION("Exception - missing .grm.id file")
    {
        Eigen::MatrixXd matrix(2, 2);
        matrix << 1.0, 0.5, 0.5, 1.0;

        grm_files.create_bin_only(matrix);

        REQUIRE_THROWS_AS(
            GrmLoader(grm_files.prefix()), gelex::FileNotFoundException);
    }
}

TEST_CASE(
    "GrmLoader - Constructor with size mismatch",
    "[grm_loader][construction][exception]")
{
    FileFixture files;
    GrmFileFixture grm_files(files);

    SECTION("Exception - bin file size doesn't match id count")
    {
        // Create a 3x3 matrix but only 2 IDs
        Eigen::MatrixXd matrix(3, 3);
        matrix << 1.0, 0.5, 0.3, 0.5, 1.0, 0.4, 0.3, 0.4, 1.0;

        std::vector<std::string> ids = {"FAM1_IND1", "FAM2_IND2"};

        grm_files.create(matrix, ids);

        REQUIRE_THROWS_AS(
            GrmLoader(grm_files.prefix()), gelex::FileFormatException);
    }

    SECTION("Exception - error message contains expected and actual size")
    {
        Eigen::MatrixXd matrix(4, 4);
        matrix.setIdentity();

        std::vector<std::string> ids = {"FAM1_IND1", "FAM2_IND2"};

        grm_files.create(matrix, ids);

        REQUIRE_THROWS_MATCHES(
            GrmLoader(grm_files.prefix()),
            gelex::FileFormatException,
            MessageMatches(ContainsSubstring("size mismatch")));
    }
}

// ============================================================================
// Accessor tests
// ============================================================================

TEST_CASE("GrmLoader - sample_ids accessor", "[grm_loader][accessor]")
{
    FileFixture files;
    GrmFileFixture grm_files(files);

    SECTION("Happy path - verify loaded sample IDs")
    {
        const Eigen::Index n = 3;
        auto matrix = make_symmetric_matrix(n);
        std::vector<std::string> ids = {"FAM1_IND1", "FAM2_IND2", "FAM3_IND3"};

        grm_files.create(matrix, ids);
        GrmLoader loader(grm_files.prefix());

        auto loaded_ids = loader.sample_ids();
        REQUIRE(loaded_ids.size() == 3);
        REQUIRE(loaded_ids[0] == "FAM1_IND1");
        REQUIRE(loaded_ids[1] == "FAM2_IND2");
        REQUIRE(loaded_ids[2] == "FAM3_IND3");
    }

    SECTION("Happy path - IDs with multiple underscores preserved")
    {
        const Eigen::Index n = 2;
        auto matrix = make_symmetric_matrix(n);
        // Writer splits "A_B_C" -> FID="A", IID="B_C"
        // Loader reads "A\tB_C" -> "A_B_C"
        std::vector<std::string> ids = {"FAM_1_IND_1", "FAM_2_IND_2"};

        grm_files.create(matrix, ids);
        GrmLoader loader(grm_files.prefix());

        auto loaded_ids = loader.sample_ids();
        REQUIRE(loaded_ids.size() == 2);
        // Writer: "FAM_1_IND_1" -> "FAM\t1_IND_1"
        // Loader: "FAM\t1_IND_1" -> "FAM_1_IND_1"
        REQUIRE(loaded_ids[0] == "FAM_1_IND_1");
        REQUIRE(loaded_ids[1] == "FAM_2_IND_2");
    }
}

TEST_CASE("GrmLoader - num_samples accessor", "[grm_loader][accessor]")
{
    FileFixture files;
    GrmFileFixture grm_files(files);

    SECTION("Happy path - verify num_samples for small matrix")
    {
        const Eigen::Index n = 5;
        auto matrix = make_symmetric_matrix(n);
        auto ids = make_sample_ids(n);

        grm_files.create(matrix, ids);
        GrmLoader loader(grm_files.prefix());

        REQUIRE(loader.num_samples() == n);
    }

    SECTION("Happy path - verify num_samples for larger matrix")
    {
        const Eigen::Index n = 50;
        auto matrix = make_symmetric_matrix(n);
        auto ids = make_sample_ids(n);

        grm_files.create(matrix, ids);
        GrmLoader loader(grm_files.prefix());

        REQUIRE(loader.num_samples() == n);
    }
}

// ============================================================================
// load() - Complete matrix loading tests
// ============================================================================

TEST_CASE("GrmLoader - Load complete 3x3 GRM", "[grm_loader][load]")
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

        GrmLoader loader(grm_files.prefix());
        auto loaded = loader.load();

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

TEST_CASE("GrmLoader - Load complete 10x10 GRM", "[grm_loader][load]")
{
    FileFixture files;
    GrmFileFixture grm_files(files);

    SECTION("Happy path - load and verify 10x10 matrix")
    {
        const Eigen::Index n = 10;
        auto original = make_symmetric_matrix(n);
        auto ids = make_sample_ids(n);

        grm_files.create(original, ids);

        GrmLoader loader(grm_files.prefix());
        auto loaded = loader.load();

        REQUIRE(loaded.rows() == n);
        REQUIRE(loaded.cols() == n);

        // Verify dimensions and some spot checks
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

TEST_CASE("GrmLoader - Verify loaded matrix is symmetric", "[grm_loader][load]")
{
    FileFixture files;
    GrmFileFixture grm_files(files);

    SECTION("Happy path - loaded matrix should be symmetric")
    {
        const Eigen::Index n = 5;
        auto original = make_symmetric_matrix(n);
        auto ids = make_sample_ids(n);

        grm_files.create(original, ids);

        GrmLoader loader(grm_files.prefix());
        auto loaded = loader.load();

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
    "GrmLoader - Verify numerical precision (double-float-double)",
    "[grm_loader][load][numerical]")
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

        GrmLoader loader(grm_files.prefix());
        auto loaded = loader.load();

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

TEST_CASE("GrmLoader - Load with subset of IDs", "[grm_loader][load_filtered]")
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
            = {"FAM0_IND0", "FAM1_IND1", "FAM2_IND2", "FAM3_IND3"};

        grm_files.create(original, ids);

        GrmLoader loader(grm_files.prefix());

        // Load only samples 1 and 3, mapping to indices 0 and 1
        std::unordered_map<std::string, Eigen::Index> id_map
            = {{"FAM1_IND1", 0}, {"FAM3_IND3", 1}};

        auto loaded = loader.load(id_map);

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

TEST_CASE("GrmLoader - Load with reordered IDs", "[grm_loader][load_filtered]")
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

        std::vector<std::string> ids = {"FAM0_IND0", "FAM1_IND1", "FAM2_IND2"};

        grm_files.create(original, ids);

        GrmLoader loader(grm_files.prefix());

        // Reverse the order: original[2]->0, original[1]->1, original[0]->2
        std::unordered_map<std::string, Eigen::Index> id_map
            = {{"FAM2_IND2", 0}, {"FAM1_IND1", 1}, {"FAM0_IND0", 2}};

        auto loaded = loader.load(id_map);

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

TEST_CASE("GrmLoader - Load with empty id_map", "[grm_loader][load_filtered]")
{
    FileFixture files;
    GrmFileFixture grm_files(files);

    SECTION("Happy path - empty id_map returns empty matrix")
    {
        const Eigen::Index n = 3;
        auto original = make_symmetric_matrix(n);
        auto ids = make_sample_ids(n);

        grm_files.create(original, ids);

        GrmLoader loader(grm_files.prefix());

        std::unordered_map<std::string, Eigen::Index> empty_map;
        auto loaded = loader.load(empty_map);

        REQUIRE(loaded.rows() == 0);
        REQUIRE(loaded.cols() == 0);
    }
}

TEST_CASE(
    "GrmLoader - Load with invalid ID throws exception",
    "[grm_loader][load_filtered][exception]")
{
    FileFixture files;
    GrmFileFixture grm_files(files);

    SECTION("Exception - ID not found in file")
    {
        const Eigen::Index n = 3;
        auto original = make_symmetric_matrix(n);
        auto ids = make_sample_ids(n);

        grm_files.create(original, ids);

        GrmLoader loader(grm_files.prefix());

        std::unordered_map<std::string, Eigen::Index> id_map
            = {{"FAM0_IND0", 0}, {"NONEXISTENT_ID", 1}};

        REQUIRE_THROWS_AS(loader.load(id_map), gelex::InvalidInputException);
    }

    SECTION("Exception - error message contains the invalid ID")
    {
        const Eigen::Index n = 2;
        auto original = make_symmetric_matrix(n);
        auto ids = make_sample_ids(n);

        grm_files.create(original, ids);

        GrmLoader loader(grm_files.prefix());

        std::unordered_map<std::string, Eigen::Index> id_map
            = {{"MISSING_SAMPLE", 0}};

        REQUIRE_THROWS_MATCHES(
            loader.load(id_map),
            gelex::InvalidInputException,
            MessageMatches(ContainsSubstring("MISSING_SAMPLE")));
    }
}

TEST_CASE(
    "GrmLoader - Load with non-contiguous target indices",
    "[grm_loader][load_filtered]")
{
    FileFixture files;
    GrmFileFixture grm_files(files);

    SECTION("Happy path - sparse target indices creates larger matrix")
    {
        Eigen::MatrixXd original(3, 3);
        // clang-format off
        original << 1.0, 0.1, 0.2,
                    0.1, 2.0, 0.3,
                    0.2, 0.3, 3.0;
        // clang-format on

        std::vector<std::string> ids = {"FAM0_IND0", "FAM1_IND1", "FAM2_IND2"};

        grm_files.create(original, ids);

        GrmLoader loader(grm_files.prefix());

        // Map to non-contiguous indices: 0 -> 0, 2 -> 5
        // Output matrix size should be max_idx + 1 = 6
        std::unordered_map<std::string, Eigen::Index> id_map
            = {{"FAM0_IND0", 0}, {"FAM2_IND2", 5}};

        auto loaded = loader.load(id_map);

        REQUIRE(loaded.rows() == 6);
        REQUIRE(loaded.cols() == 6);

        auto to_float = [](double v)
        { return static_cast<double>(static_cast<float>(v)); };

        // Check mapped values
        REQUIRE(loaded(0, 0) == to_float(1.0));  // original(0,0)
        REQUIRE(loaded(5, 5) == to_float(3.0));  // original(2,2)
        REQUIRE(loaded(0, 5) == to_float(0.2));  // original(0,2)
        REQUIRE(loaded(5, 0) == to_float(0.2));  // original(2,0)

        // Unmapped indices should be zero
        REQUIRE(loaded(1, 1) == 0.0);
        REQUIRE(loaded(2, 2) == 0.0);
        REQUIRE(loaded(3, 3) == 0.0);
        REQUIRE(loaded(4, 4) == 0.0);
    }
}

// ============================================================================
// ID parsing tests
// ============================================================================

TEST_CASE("GrmLoader - ID parsing from file", "[grm_loader][id_parsing]")
{
    FileFixture files;
    GrmFileFixture grm_files(files);

    SECTION("Happy path - IDs without underscore in original become duplicated")
    {
        // When GrmIdWriter writes ID without underscore, it writes "ID\tID"
        // GrmLoader reads "ID\tID" as "ID_ID"
        const Eigen::Index n = 2;
        auto matrix = make_symmetric_matrix(n);
        std::vector<std::string> ids = {"SAMPLE1", "SAMPLE2"};

        grm_files.create(matrix, ids);
        GrmLoader loader(grm_files.prefix());

        auto loaded_ids = loader.sample_ids();
        // Writer: "SAMPLE1" -> "SAMPLE1\tSAMPLE1"
        // Loader: "SAMPLE1\tSAMPLE1" -> "SAMPLE1_SAMPLE1"
        REQUIRE(loaded_ids[0] == "SAMPLE1_SAMPLE1");
        REQUIRE(loaded_ids[1] == "SAMPLE2_SAMPLE2");
    }
}

// ============================================================================
// Round-trip verification tests
// ============================================================================

TEST_CASE(
    "GrmLoader - Round-trip write and load verification",
    "[grm_loader][roundtrip]")
{
    FileFixture files;
    GrmFileFixture grm_files(files);

    SECTION("Happy path - verify data integrity through write/load cycle")
    {
        const Eigen::Index n = 20;
        auto original = make_symmetric_matrix(n);
        auto ids = make_sample_ids(n);

        grm_files.create(original, ids);

        GrmLoader loader(grm_files.prefix());
        auto loaded = loader.load();

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

        GrmLoader loader(grm_files.prefix());
        auto loaded = loader.load();

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
    "GrmLoader - Load matrix with special values",
    "[grm_loader][special]")
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

        GrmLoader loader(grm_files.prefix());
        auto loaded = loader.load();

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

        GrmLoader loader(grm_files.prefix());
        auto loaded = loader.load();

        REQUIRE(std::isnan(loaded(0, 0)));
    }
}
