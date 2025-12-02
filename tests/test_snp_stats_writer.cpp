#include <algorithm>
#include <array>
#include <cstdint>
#include <fstream>
#include <vector>

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_exception.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>

#include "file_fixture.h"
#include "gelex/exception.h"

#include "../src/data/snp_stats_writer.h"

namespace fs = std::filesystem;

using namespace gelex::detail;  // NOLINT
using Catch::Matchers::EndsWith;
using gelex::test::FileFixture;

TEST_CASE("SnpStatsWriter - Constructor and path access", "[data][snp_stats]")
{
    FileFixture files;

    SECTION("Happy path - create writer with valid path")
    {
        auto file_path = files.create_empty_file(".snp_stats");

        REQUIRE_NOTHROW(
            [&]()
            {
                SnpStatsWriter writer(file_path);
                REQUIRE(writer.path() == file_path);
            }());
    }
}

TEST_CASE("SnpStatsWriter - Write valid data", "[data][snp_stats]")
{
    FileFixture files;

    SECTION("Happy path - write basic data with monomorphic variants")
    {
        auto file_path = files.create_empty_file(".snp_stats");

        const int64_t num_samples = 100;
        const std::vector<int64_t> monomorphic_indices = {2, 5, 8};
        const std::vector<double> means
            = {0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0};
        const std::vector<double> stddevs
            = {0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5};

        REQUIRE_NOTHROW(
            [&]()
            {
                SnpStatsWriter writer(file_path);
                writer.write(num_samples, monomorphic_indices, means, stddevs);
            }());

        std::ifstream file(file_path, std::ios::binary);
        REQUIRE(file.is_open());

        std::array<int64_t, 3> header{};
        file.read(reinterpret_cast<char*>(header.data()), sizeof(int64_t) * 3);
        REQUIRE(header[0] == num_samples);
        REQUIRE(header[1] == static_cast<int64_t>(means.size()));
        REQUIRE(header[2] == static_cast<int64_t>(monomorphic_indices.size()));

        std::vector<int64_t> read_monomorphic(monomorphic_indices.size());
        if (!monomorphic_indices.empty())
        {
            file.read(
                reinterpret_cast<char*>(read_monomorphic.data()),
                static_cast<std::streamsize>(
                    monomorphic_indices.size() * sizeof(int64_t)));
            REQUIRE(read_monomorphic == monomorphic_indices);
        }

        std::vector<double> read_means(means.size());
        file.read(
            reinterpret_cast<char*>(read_means.data()),
            static_cast<std::streamsize>(means.size() * sizeof(double)));
        REQUIRE(std::equal(means.begin(), means.end(), read_means.begin()));

        std::vector<double> read_stddevs(stddevs.size());
        file.read(
            reinterpret_cast<char*>(read_stddevs.data()),
            static_cast<std::streamsize>(stddevs.size() * sizeof(double)));
        REQUIRE(
            std::equal(stddevs.begin(), stddevs.end(), read_stddevs.begin()));
    }

    SECTION("Happy path - write data without monomorphic variants")
    {
        auto file_path = files.create_empty_file(".snp_stats");

        const int64_t num_samples = 50;
        const std::vector<int64_t> monomorphic_indices = {};
        const std::vector<double> means = {0.5, 0.6, 0.7};
        const std::vector<double> stddevs = {0.1, 0.2, 0.3};

        REQUIRE_NOTHROW(
            [&]()
            {
                SnpStatsWriter writer(file_path);
                writer.write(num_samples, monomorphic_indices, means, stddevs);
            }());

        std::ifstream file(file_path, std::ios::binary);
        REQUIRE(file.is_open());

        std::array<int64_t, 3> header{};
        file.read(reinterpret_cast<char*>(header.data()), sizeof(int64_t) * 3);
        REQUIRE(header[0] == num_samples);
        REQUIRE(header[1] == 3);
        REQUIRE(header[2] == 0);

        std::vector<double> read_means(3);
        file.read(
            reinterpret_cast<char*>(read_means.data()), 3 * sizeof(double));
        REQUIRE(read_means == means);
    }

    SECTION("Happy path - write single variant")
    {
        auto file_path = files.create_empty_file(".snp_stats");

        const int64_t num_samples = 10;
        const std::vector<int64_t> monomorphic_indices = {0};
        const std::vector<double> means = {0.8};
        const std::vector<double> stddevs = {0.2};

        REQUIRE_NOTHROW(
            [&]()
            {
                SnpStatsWriter writer(file_path);
                writer.write(num_samples, monomorphic_indices, means, stddevs);
            }());

        std::ifstream file(file_path, std::ios::binary);
        std::array<int64_t, 3> header{};
        file.read(reinterpret_cast<char*>(header.data()), sizeof(int64_t) * 3);
        REQUIRE(header[0] == num_samples);
        REQUIRE(header[1] == 1);
        REQUIRE(header[2] == 1);
    }
}

TEST_CASE("SnpStatsWriter - Argument validation", "[data][snp_stats]")
{
    FileFixture files;

    SECTION("Exception - means and stddevs size mismatch")
    {
        auto file_path = files.create_empty_file(".snp_stats");

        const int64_t num_samples = 100;
        const std::vector<int64_t> monomorphic_indices = {};
        const std::vector<double> means = {0.1, 0.2, 0.3};
        const std::vector<double> stddevs = {0.1, 0.2};

        SnpStatsWriter writer(file_path);
        REQUIRE_THROWS_MATCHES(
            writer.write(num_samples, monomorphic_indices, means, stddevs),
            gelex::ArgumentValidationException,
            Catch::Matchers::MessageMatches(EndsWith(
                "means (3) and stddevs (2) must have the same length.")));
    }

    SECTION("Exception - empty means array")
    {
        auto file_path = files.create_empty_file(".snp_stats");

        const int64_t num_samples = 100;
        const std::vector<int64_t> monomorphic_indices = {};
        const std::vector<double> means = {};
        const std::vector<double> stddevs = {};

        SnpStatsWriter writer(file_path);
        REQUIRE_THROWS_MATCHES(
            writer.write(num_samples, monomorphic_indices, means, stddevs),
            gelex::ArgumentValidationException,
            Catch::Matchers::MessageMatches(
                EndsWith("means and stddevs cannot be empty")));
    }

    SECTION("Exception - empty stddevs array")
    {
        auto file_path = files.create_empty_file(".snp_stats");

        const int64_t num_samples = 100;
        const std::vector<int64_t> monomorphic_indices = {};
        const std::vector<double> means = {0.1};
        const std::vector<double> stddevs = {};

        SnpStatsWriter writer(file_path);
        REQUIRE_THROWS_MATCHES(
            writer.write(num_samples, monomorphic_indices, means, stddevs),
            gelex::ArgumentValidationException,
            Catch::Matchers::MessageMatches(EndsWith(
                "means (1) and stddevs (0) must have the same length.")));
    }

    SECTION("Exception - monomorphic index out of range")
    {
        auto file_path = files.create_empty_file(".snp_stats");

        const int64_t num_samples = 100;
        const std::vector<int64_t> monomorphic_indices = {5};
        const std::vector<double> means = {0.1, 0.2, 0.3};
        const std::vector<double> stddevs = {0.1, 0.2, 0.3};

        SnpStatsWriter writer(file_path);
        REQUIRE_THROWS_MATCHES(
            writer.write(num_samples, monomorphic_indices, means, stddevs),
            gelex::ArgumentValidationException,
            Catch::Matchers::MessageMatches(
                EndsWith("Monomorphic SNP index 5 is out of range [0, 3).")));
    }
}
