#include <filesystem>
#include <fstream>
#include <iostream>
#include <vector>

#include <catch2/catch_test_macros.hpp>

#include "../src/data/snp_stats_writer.h"
#include "gelex/error.h"

class SnpStatsWriterTestFixture
{
   public:
    SnpStatsWriterTestFixture()
    {
        createValidTestFile();
        createReadOnlyFile();
    }

    SnpStatsWriterTestFixture(const SnpStatsWriterTestFixture&) = default;
    SnpStatsWriterTestFixture(SnpStatsWriterTestFixture&&) = delete;
    SnpStatsWriterTestFixture& operator=(const SnpStatsWriterTestFixture&)
        = default;
    SnpStatsWriterTestFixture& operator=(SnpStatsWriterTestFixture&&) = delete;
    ~SnpStatsWriterTestFixture() { remove(); }

    static void remove()
    {
        std::remove("test_valid.snpstats");
        std::remove("test_readonly.snpstats");
        std::remove("test_output.snpstats");
    }

    static void createValidTestFile()
    {
        // Create an empty valid file for testing
        std::ofstream file("test_valid.snpstats");
    }

    static void createReadOnlyFile()
    {
        // Create a file and make it read-only
        std::ofstream file("test_readonly.snpstats");
        file.close();
        std::filesystem::permissions(
            "test_readonly.snpstats", std::filesystem::perms::owner_read);
    }

    static std::vector<int64_t> createMonomorphicIndices() { return {2, 5, 8}; }

    static std::vector<double> createMeans()
    {
        return {0.1, 0.2, 0.3, 0.4, 0.5};
    }

    static std::vector<double> createStddevs()
    {
        return {0.05, 0.06, 0.07, 0.08, 0.09};
    }
};

TEST_CASE_PERSISTENT_FIXTURE(
    SnpStatsWriterTestFixture,
    "SnpStatsWriter::create function",
    "[snp_stats][writer]")
{
    SECTION("Valid file path")
    {
        auto result
            = gelex::detail::SnpStatsWriter::create("test_valid.snpstats");
        REQUIRE(result.has_value());
    }

    SECTION("Non-existent file")
    {
        auto result = gelex::detail::SnpStatsWriter::create(
            "non_existent_file.snpstats");
        REQUIRE(result.has_value());  // Should create the file
    }

    SECTION("Read-only file")
    {
        auto result
            = gelex::detail::SnpStatsWriter::create("test_readonly.snpstats");
        REQUIRE_FALSE(result.has_value());
        REQUIRE(result.error().code == gelex::ErrorCode::FileIOError);
    }
}

TEST_CASE_PERSISTENT_FIXTURE(
    SnpStatsWriterTestFixture,
    "SnpStatsWriter::write function - valid data",
    "[snp_stats][writer]")
{
    SECTION("Write valid data")
    {
        auto writer_result
            = gelex::detail::SnpStatsWriter::create("test_output.snpstats");
        REQUIRE(writer_result.has_value());

        auto& writer = writer_result.value();

        auto monomorphic_indices
            = SnpStatsWriterTestFixture::createMonomorphicIndices();
        auto means = SnpStatsWriterTestFixture::createMeans();
        auto stddevs = SnpStatsWriterTestFixture::createStddevs();

        auto write_result = writer.write(
            100,  // num_samples
            5,    // num_variants
            3,    // num_monomorphic
            monomorphic_indices,
            means,
            stddevs);

        REQUIRE(write_result.has_value());

        // Verify file was written correctly
        std::ifstream file("test_output.snpstats", std::ios::binary);
        REQUIRE(file.is_open());

        // Read header
        int64_t num_samples = 0;
        int64_t num_variants = 0;
        int64_t num_monomorphic = 0;
        file.read(reinterpret_cast<char*>(&num_samples), sizeof(int64_t));
        file.read(reinterpret_cast<char*>(&num_variants), sizeof(int64_t));
        file.read(reinterpret_cast<char*>(&num_monomorphic), sizeof(int64_t));

        REQUIRE(num_samples == 100);
        REQUIRE(num_variants == 5);
        REQUIRE(num_monomorphic == 3);

        // Read monomorphic indices
        std::vector<int64_t> read_indices(num_monomorphic);
        file.read(
            reinterpret_cast<char*>(read_indices.data()),
            static_cast<std::streamsize>(num_monomorphic * sizeof(int64_t)));
        REQUIRE(read_indices == monomorphic_indices);

        // Read means
        std::vector<double> read_means(num_variants);
        file.read(
            reinterpret_cast<char*>(read_means.data()),
            static_cast<std::streamsize>(num_variants * sizeof(double)));
        REQUIRE(read_means == means);

        // Read stddevs
        std::vector<double> read_stddevs(num_variants);
        file.read(
            reinterpret_cast<char*>(read_stddevs.data()),
            static_cast<std::streamsize>(num_variants * sizeof(double)));
        REQUIRE(read_stddevs == stddevs);
    }
}

TEST_CASE_PERSISTENT_FIXTURE(
    SnpStatsWriterTestFixture,
    "SnpStatsWriter::write function - error conditions",
    "[snp_stats][writer]")
{
    SECTION("Size mismatch - num_monomorphic vs indices")
    {
        auto writer_result
            = gelex::detail::SnpStatsWriter::create("test_output.snpstats");
        REQUIRE(writer_result.has_value());

        auto& writer = writer_result.value();

        auto write_result = writer.write(
            100,        // num_samples
            5,          // num_variants
            2,          // num_monomorphic (should be 3)
            {2, 5, 8},  // 3 indices
            {0.1, 0.2, 0.3, 0.4, 0.5},
            {0.05, 0.06, 0.07, 0.08, 0.09});

        REQUIRE_FALSE(write_result.has_value());
        REQUIRE(write_result.error().code == gelex::ErrorCode::InvalidArgument);
    }

    SECTION("Size mismatch - means/stddevs vs num_variants")
    {
        auto writer_result
            = gelex::detail::SnpStatsWriter::create("test_output.snpstats");
        REQUIRE(writer_result.has_value());

        auto& writer = writer_result.value();

        auto write_result = writer.write(
            100,  // num_samples
            5,    // num_variants
            3,    // num_monomorphic
            {2, 5, 8},
            {0.1, 0.2, 0.3},  // Only 3 means (should be 5)
            {0.05, 0.06, 0.07, 0.08, 0.09});

        REQUIRE_FALSE(write_result.has_value());
        REQUIRE(write_result.error().code == gelex::ErrorCode::InvalidArgument);
    }

    SECTION("Empty monomorphic indices should be allowed")
    {
        auto writer_result
            = gelex::detail::SnpStatsWriter::create("test_output.snpstats");
        REQUIRE(writer_result.has_value());

        auto& writer = writer_result.value();

        auto write_result = writer.write(
            100,  // num_samples
            5,    // num_variants
            0,    // num_monomorphic
            {},   // empty indices (should be allowed)
            {0.1, 0.2, 0.3, 0.4, 0.5},
            {0.05, 0.06, 0.07, 0.08, 0.09});

        REQUIRE(write_result.has_value());
    }

    SECTION("All empty vectors")
    {
        auto writer_result
            = gelex::detail::SnpStatsWriter::create("test_output.snpstats");
        REQUIRE(writer_result.has_value());

        auto& writer = writer_result.value();

        auto write_result = writer.write(
            100,  // num_samples
            0,    // num_variants
            0,    // num_monomorphic
            {},   // empty indices
            {},   // empty means
            {}    // empty stddevs
        );

        REQUIRE_FALSE(write_result.has_value());
        REQUIRE(write_result.error().code == gelex::ErrorCode::InvalidArgument);
    }
}

TEST_CASE_PERSISTENT_FIXTURE(
    SnpStatsWriterTestFixture,
    "SnpStatsWriter file format verification",
    "[snp_stats][writer]")
{
    SECTION("Verify binary format structure")
    {
        auto writer_result
            = gelex::detail::SnpStatsWriter::create("test_output.snpstats");
        REQUIRE(writer_result.has_value());

        auto& writer = writer_result.value();

        std::vector<int64_t> monomorphic_indices = {1, 3, 7};
        std::vector<double> means = {0.15, 0.25, 0.35, 0.45};
        std::vector<double> stddevs = {0.055, 0.065, 0.075, 0.085};

        auto write_result = writer.write(
            50,  // num_samples
            4,   // num_variants
            3,   // num_monomorphic
            monomorphic_indices,
            means,
            stddevs);

        REQUIRE(write_result.has_value());

        // Verify file size is correct
        std::ifstream file(
            "test_output.snpstats", std::ios::binary | std::ios::ate);
        auto file_size = file.tellg();

        size_t expected_size = (3 * sizeof(int64_t)) +  // header
                               (3 * sizeof(int64_t)) +  // monomorphic indices
                               (4 * sizeof(double)) +   // means
                               (4 * sizeof(double));    // stddevs

        REQUIRE(static_cast<size_t>(file_size) == expected_size);
    }
}
