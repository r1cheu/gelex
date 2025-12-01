#include <catch2/catch_test_macros.hpp>

#include "file_fixture.h"
#include "gelex/exception.h"

#include "../src/data/loader/bim_loader.h"

using namespace gelex::detail;  // NOLINT
using gelex::test::FileFixture;

TEST_CASE("BimLoader - Valid file parsing", "[loader][bim]")
{
    FileFixture files;

    SECTION("Happy path - tab-delimited file")
    {
        auto file_path = files.create_text_file(
            "1\trs12345\t0\t1000\tA\tG\n"
            "1\trs67890\t0.001\t2000\tC\tT\n"
            "2\trs24680\t0.002\t3000\tG\tA\n"
            "X\trs13579\t0.003\t4000\tT\tC\n"
            "1\trs11223\t0.004\t5000\tA\tT",
            ".bim");

        REQUIRE_NOTHROW(
            [&]()
            {
                BimLoader bim_loader(file_path);

                const auto snp_ids = bim_loader.get_ids();
                REQUIRE(snp_ids.size() == 5);
                REQUIRE(snp_ids[0] == "rs12345");
                REQUIRE(snp_ids[1] == "rs67890");
                REQUIRE(snp_ids[2] == "rs24680");
                REQUIRE(snp_ids[3] == "rs13579");
                REQUIRE(snp_ids[4] == "rs11223");

                const auto& snp1 = bim_loader[0];
                REQUIRE(snp1.chrom == '1');
                REQUIRE(snp1.id == "rs12345");
                REQUIRE(snp1.base_coordinate == 1000);
                REQUIRE(snp1.A1 == 'A');
                REQUIRE(snp1.A2 == 'G');

                const auto& snp3 = bim_loader[3];
                REQUIRE(snp3.chrom == 'X');
                REQUIRE(snp3.id == "rs13579");
                REQUIRE(snp3.base_coordinate == 4000);
                REQUIRE(snp3.A1 == 'T');
                REQUIRE(snp3.A2 == 'C');
            }());
    }

    SECTION("Happy path - space-delimited file")
    {
        auto file_path = files.create_text_file(
            "1 rs12345 0 1000 A G\n"
            "1 rs67890 0.001 2000 C T\n"
            "2 rs24680 0.002 3000 G A",
            ".bim");

        REQUIRE_NOTHROW(
            [&]()
            {
                BimLoader bim_loader(file_path);

                const auto& snp_ids = bim_loader.get_ids();
                REQUIRE(snp_ids.size() == 3);
                REQUIRE(snp_ids[0] == "rs12345");
                REQUIRE(snp_ids[1] == "rs67890");
                REQUIRE(snp_ids[2] == "rs24680");

                const auto& snp1 = bim_loader[0];
                REQUIRE(snp1.chrom == '1');
                REQUIRE(snp1.id == "rs12345");
                REQUIRE(snp1.base_coordinate == 1000);
                REQUIRE(snp1.A1 == 'A');
                REQUIRE(snp1.A2 == 'G');

                const auto& snp3 = bim_loader[2];
                REQUIRE(snp3.chrom == '2');
                REQUIRE(snp3.id == "rs24680");
                REQUIRE(snp3.base_coordinate == 3000);
                REQUIRE(snp3.A1 == 'G');
                REQUIRE(snp3.A2 == 'A');
            }());
    }
}

TEST_CASE("BimLoader - Malformed column count", "[loader][bim]")
{
    FileFixture files;

    SECTION("Exception - inconsistent columns")
    {
        auto file_path = files.create_text_file(
            "1\trs12345\t0\t1000\tA\tG\n"
            "1\trs67890\t0.001\t2000\tC",  // Missing last column
            ".bim");
        REQUIRE_THROWS_AS(BimLoader(file_path), gelex::FileFormatException);
    }
}

TEST_CASE("BimLoader - Insufficient columns", "[loader][bim]")
{
    FileFixture files;

    SECTION("Exception - single column")
    {
        auto file_path = files.create_text_file("1", ".bim");
        REQUIRE_THROWS_AS(BimLoader(file_path), gelex::FileFormatException);
    }

    SECTION("Exception - five columns")
    {
        auto file_path
            = files.create_text_file("1\trs12345\t0\t1000\tA", ".bim");
        REQUIRE_THROWS_AS(BimLoader(file_path), gelex::FileFormatException);
    }
}

TEST_CASE("BimLoader - Invalid position data", "[loader][bim]")
{
    FileFixture files;

    SECTION("Exception - non-numeric position")
    {
        auto file_path
            = files.create_text_file("1\trs12345\t0\tinvalid\tA\tG", ".bim");
        REQUIRE_THROWS_AS(BimLoader(file_path), gelex::FileFormatException);
    }

    SECTION("Exception - empty position field")
    {
        auto file_path
            = files.create_text_file("1\trs12345\t0\t\tA\tG", ".bim");
        REQUIRE_THROWS_AS(BimLoader(file_path), gelex::FileFormatException);
    }
}

TEST_CASE("BimLoader - Comprehensive happy path tests", "[loader][bim]")
{
    FileFixture files;

    SECTION("Test all public methods")
    {
        auto file_path = files.create_text_file(
            "1\trs12345\t0\t1000\tA\tG\n"
            "2\trs67890\t0.001\t2000\tC\tT",
            ".bim");

        BimLoader bim_loader(file_path);

        // Test meta() method
        const auto& meta = bim_loader.info();
        REQUIRE(meta.size() == 2);
        REQUIRE(meta[0].id == "rs12345");
        REQUIRE(meta[1].id == "rs67890");

        // Test get_ids() method
        const auto ids = bim_loader.get_ids();
        REQUIRE(ids.size() == 2);
        REQUIRE(ids[0] == "rs12345");
        REQUIRE(ids[1] == "rs67890");

        // Test operator[]
        REQUIRE_NOTHROW(bim_loader[0]);
        REQUIRE(bim_loader[0].id == "rs12345");

        // Test take_meta() method
        auto moved_meta = std::move(bim_loader).take_info();
        REQUIRE(moved_meta.size() == 2);
        REQUIRE(moved_meta[0].id == "rs12345");
        REQUIRE(moved_meta[1].id == "rs67890");
    }
}
