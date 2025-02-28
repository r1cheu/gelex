#include <string>
#include <vector>

#include <armadillo>
#include <catch2/catch_test_macros.hpp>
#include "chenx/data/bed_reader.h"

TEST_CASE("BedReader initialization", "[bedreader]")
{
    const std::string test_bed
        = std::string(CHENX_TESTS_DIR) + "/data/test.bed";
    const size_t chunk_size = 10;

    SECTION("Valid file initialization")
    {
        REQUIRE_NOTHROW(chenx::BedReader(test_bed, {}, chunk_size));

        chenx::BedReader reader(test_bed, {}, chunk_size);
        REQUIRE(reader.num_snps() > 0);
        REQUIRE(reader.num_individuals() > 0);
    }

    SECTION("Invalid file throws")
    {
        REQUIRE_THROWS_AS(
            chenx::BedReader("invalid_path.bed", {}), std::runtime_error);
    }
}

TEST_CASE("BedReader small chunk reading", "[bedreader]")
{
    const std::string test_bed
        = std::string(CHENX_TESTS_DIR) + "/data/test.bed";
    chenx::BedReader reader(test_bed, {}, 2);

    SECTION("HasNext returns correct state")
    {
        REQUIRE(reader.HasNext());

        // Read all chunks
        while (reader.HasNext())
        {
            reader.ReadChunk();
        }

        REQUIRE_FALSE(reader.HasNext());
    }

    SECTION("ReadChunk returns correct matrix dimensions")
    {
        auto chunk = reader.ReadChunk();
        REQUIRE(chunk.n_rows == reader.num_individuals());
        REQUIRE(chunk.n_cols <= 10);  // chunk_size
    }
}

TEST_CASE("BedReader Big chunk reading", "[bedreader]")
{
    const std::string test_bed
        = std::string(CHENX_TESTS_DIR) + "/data/test.bed";
    chenx::BedReader reader(test_bed, {}, 10);

    SECTION("HasNext returns correct state")
    {
        REQUIRE(reader.HasNext());

        // Read all chunks
        while (reader.HasNext())
        {
            reader.ReadChunk();
        }

        REQUIRE_FALSE(reader.HasNext());
    }

    SECTION("ReadChunk returns correct matrix dimensions")
    {
        auto chunk = reader.ReadChunk();
        REQUIRE(chunk.n_rows == reader.num_individuals());
        REQUIRE(chunk.n_cols <= 10);  // chunk_size
    }
}

TEST_CASE("BedReader metadata access", "[bedreader]")
{
    const std::string test_bed
        = std::string(CHENX_TESTS_DIR) + "/data/test.bed";
    chenx::BedReader reader(test_bed, {}, 2);

    SECTION("SNP access")
    {
        auto snps = reader.snps();
        REQUIRE_FALSE(snps.empty());
        REQUIRE(snps.size() == reader.num_snps());
    }

    SECTION("Individual access")
    {
        auto individuals = reader.individuals();
        REQUIRE_FALSE(individuals.empty());
        REQUIRE(individuals.size() == reader.num_individuals());
    }
}

TEST_CASE("BedReader exclude individuals", "[bedreader]")
{
    const std::string test_bed
        = std::string(CHENX_TESTS_DIR) + "/data/test.bed";
    std::vector<std::string> dropped_individuals{"iid2"};
    chenx::BedReader reader(test_bed, dropped_individuals, 2);

    REQUIRE(reader.num_individuals() == 2);
    for (const auto& individual : reader.individuals())
    {
        REQUIRE(individual != "iid2");
    }
    while (reader.HasNext())
    {
        arma::dmat chunk{reader.ReadChunk()};
        REQUIRE(chunk.n_rows == 2);
    }

    chenx::BedReader reader2(test_bed, dropped_individuals, 10);
    arma::dmat chunk{reader2.ReadChunk()};
    REQUIRE(chunk.n_rows == 2);
}
