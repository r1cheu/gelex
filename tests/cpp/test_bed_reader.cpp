#include <string>
#include <vector>

#include <Eigen/Core>
#include <catch2/catch_test_macros.hpp>

#include "gelex/data/bed_reader.h"

namespace test
{
}
static constexpr size_t SMALL_CHUNK_SIZE = 2;
static constexpr size_t CHUNK_SIZE = 10;

TEST_CASE("BedReader initialization", "[bedreader]")
{
    const std::string test_bed
        = std::string(GELEX_TESTS_DIR) + "/data/train.bed";

    SECTION("Valid file initialization")
    {
        REQUIRE_NOTHROW(gelex::BedReader(test_bed, CHUNK_SIZE));

        gelex::BedReader reader(test_bed, CHUNK_SIZE);
        REQUIRE(reader.num_snps() > 0);
        REQUIRE(reader.num_individuals() > 0);
    }

    SECTION("Invalid file throws")
    {
        REQUIRE_THROWS_AS(
            gelex::BedReader("invalid_path.bed"), std::runtime_error);
    }
}

TEST_CASE("BedReader small chunk reading", "[bedreader]")
{
    const std::string test_bed
        = std::string(GELEX_TESTS_DIR) + "/data/train.bed";
    gelex::BedReader reader(test_bed, SMALL_CHUNK_SIZE);

    SECTION("has_next returns correct state")
    {
        REQUIRE(reader.has_next());

        // Read all chunks
        while (reader.has_next())
        {
            reader.read_chunk();
        }

        REQUIRE_FALSE(reader.has_next());
    }

    SECTION("read_chunk returns correct matrix dimensions")
    {
        auto chunk = reader.read_chunk();
        REQUIRE(chunk.rows() == reader.num_individuals());
        REQUIRE(chunk.cols() <= 2);  // chunk_size
    }
}

TEST_CASE("BedReader Big chunk reading", "[bedreader]")
{
    const std::string test_bed
        = std::string(GELEX_TESTS_DIR) + "/data/train.bed";
    gelex::BedReader reader(test_bed, CHUNK_SIZE);

    SECTION("has_next returns correct state")
    {
        REQUIRE(reader.has_next());

        // Read all chunks
        while (reader.has_next())
        {
            reader.read_chunk();
        }

        REQUIRE_FALSE(reader.has_next());
    }

    SECTION("read_chunk returns correct matrix dimensions")
    {
        auto chunk = reader.read_chunk();
        REQUIRE(chunk.rows() == reader.num_individuals());
        REQUIRE(chunk.cols() <= 10);  // chunk_size
    }
}

TEST_CASE("BedReader metadata access", "[bedreader]")
{
    const std::string test_bed
        = std::string(GELEX_TESTS_DIR) + "/data/train.bed";
    gelex::BedReader reader(test_bed, SMALL_CHUNK_SIZE);

    SECTION("SNP access")
    {
        std::vector<std::string> expect_snps{"sid1", "sid2", "sid3", "sid4"};
        auto snps = reader.snps();
        for (size_t i = 0; i < snps.size(); ++i)
        {
            REQUIRE(snps[i] == expect_snps[i]);
        }
    }

    SECTION("Individual access")
    {
        std::vector<std::string> expect_individuals{"iid1", "iid2", "iid3"};
        auto individuals = reader.individuals();
        for (size_t i = 0; i < individuals.size(); ++i)
        {
            REQUIRE(individuals[i] == expect_individuals[i]);
        }
    }
}

TEST_CASE("BedReader exclude individuals", "[bedreader]")
{
    const std::string test_bed
        = std::string(GELEX_TESTS_DIR) + "/data/test_10.bed";
    std::vector<std::string> target_order {"iid1", "iid9", "iid6", "iid7", "iid8", "iid4", "iid10"};
    gelex::BedReader reader(test_bed, SMALL_CHUNK_SIZE, target_order);

    REQUIRE(reader.num_individuals() == 7);
    for (const auto& individual : reader.individuals())
    {
        REQUIRE(individual != "iid2");
        REQUIRE(individual != "iid3");
        REQUIRE(individual != "iid5");
    }
    while (reader.has_next())
    {
        Eigen::MatrixXd chunk = reader.read_chunk();
        REQUIRE(chunk.rows() == 7);
    }

    gelex::BedReader reader2(test_bed, CHUNK_SIZE, target_order);
    Eigen::MatrixXd chunk = reader2.read_chunk();
    REQUIRE(chunk.rows() == 7);
    
    // Expected values for comparison using Eigen initializer list
    Eigen::MatrixXd expect {{
        {1.0, 0.0, 1.0, 1.0, 0.0, 1.0, 0.0, 1.0, 1.0, 0.0},
        {1.0, 0.0, 1.0, 2.0, 1.0, 2.0, 0.0, 1.0, 1.0, 2.0},
        {0.0, 1.0, 1.0, 2.0, 1.0, 2.0, 2.0, 0.0, 2.0, 1.0},
        {0.0, 1.0, 2.0, 1.0, 0.0, 1.0, 0.0, 0.0, 1.0, 2.0},
        {1.0, 2.0, 0.0, 1.0, 1.0, 1.0, 2.0, 2.0, 0.0, 1.0},
        {1.0, 1.0, 2.0, 0.0, 1.0, 0.0, 0.0, 2.0, 1.0, 0.0},
        {1.0, 1.0, 2.0, 0.0, 2.0, 0.0, 1.0, 1.0, 0.0, 1.0}
    }};
    
    // Compare matrices with tolerance using Eigen's isApprox
    REQUIRE(chunk.isApprox(expect, 1e-5));
}
