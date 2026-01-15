#include <filesystem>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_exception.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>

#include "../src/data/grm_id_writer.h"
#include "file_fixture.h"

namespace fs = std::filesystem;

using namespace gelex::detail;  // NOLINT
using Catch::Matchers::ContainsSubstring;
using Catch::Matchers::MessageMatches;
using gelex::test::FileFixture;

// Helper function to read file content as string
auto read_file_content(const fs::path& file_path) -> std::string
{
    std::ifstream file(file_path);
    std::stringstream buffer;
    buffer << file.rdbuf();
    return buffer.str();
}

// Helper function to read file lines
auto read_file_lines(const fs::path& file_path) -> std::vector<std::string>
{
    std::ifstream file(file_path);
    std::vector<std::string> lines;
    std::string line;
    while (std::getline(file, line))
    {
        lines.push_back(line);
    }
    return lines;
}

// ============================================================================
// Constructor tests
// ============================================================================

TEST_CASE(
    "GrmIdWriter - Constructor and path access",
    "[grm_id_writer][construction]")
{
    FileFixture files;

    SECTION("Happy path - create writer with valid path")
    {
        auto file_path = files.generate_random_file_path(".grm.id");

        REQUIRE_NOTHROW(
            [&]()
            {
                GrmIdWriter writer(file_path);
                REQUIRE(writer.path() == file_path);
            }());
    }

    SECTION("Happy path - file is created on construction")
    {
        auto file_path = files.generate_random_file_path(".grm.id");

        {
            GrmIdWriter writer(file_path);
        }

        REQUIRE(fs::exists(file_path));
    }
}

// ============================================================================
// Empty input tests
// ============================================================================

TEST_CASE("GrmIdWriter - Write empty ID list", "[grm_id_writer][empty]")
{
    FileFixture files;

    SECTION("Happy path - write empty vector produces empty file")
    {
        auto file_path = files.generate_random_file_path(".grm.id");

        std::vector<std::string> empty_ids;

        REQUIRE_NOTHROW(
            [&]()
            {
                GrmIdWriter writer(file_path);
                writer.write(empty_ids);
            }());

        REQUIRE(fs::exists(file_path));
        REQUIRE(fs::file_size(file_path) == 0);
    }
}

// ============================================================================
// Basic write tests
// ============================================================================

TEST_CASE("GrmIdWriter - Write single ID", "[grm_id_writer][basic]")
{
    FileFixture files;

    SECTION("Happy path - write single ID with underscore")
    {
        auto file_path = files.generate_random_file_path(".grm.id");

        std::vector<std::string> ids = {"FAM1_IND1"};

        {
            GrmIdWriter writer(file_path);
            writer.write(ids);
        }

        auto content = read_file_content(file_path);
        REQUIRE(content == "FAM1\tIND1\n");
    }
}

TEST_CASE("GrmIdWriter - Write multiple IDs", "[grm_id_writer][basic]")
{
    FileFixture files;

    SECTION("Happy path - write multiple IDs")
    {
        auto file_path = files.generate_random_file_path(".grm.id");

        std::vector<std::string> ids = {"FAM1_IND1", "FAM2_IND2", "FAM3_IND3"};

        {
            GrmIdWriter writer(file_path);
            writer.write(ids);
        }

        auto lines = read_file_lines(file_path);
        REQUIRE(lines.size() == 3);
        REQUIRE(lines[0] == "FAM1\tIND1");
        REQUIRE(lines[1] == "FAM2\tIND2");
        REQUIRE(lines[2] == "FAM3\tIND3");
    }
}

// ============================================================================
// ID splitting tests
// ============================================================================

TEST_CASE("GrmIdWriter - ID with single underscore", "[grm_id_writer][split]")
{
    FileFixture files;

    SECTION("Happy path - split by first underscore")
    {
        auto file_path = files.generate_random_file_path(".grm.id");

        std::vector<std::string> ids = {"FAMILY_INDIVIDUAL"};

        {
            GrmIdWriter writer(file_path);
            writer.write(ids);
        }

        auto content = read_file_content(file_path);
        REQUIRE(content == "FAMILY\tINDIVIDUAL\n");
    }
}

TEST_CASE(
    "GrmIdWriter - ID with multiple underscores",
    "[grm_id_writer][split]")
{
    FileFixture files;

    SECTION("Happy path - split only by first underscore")
    {
        auto file_path = files.generate_random_file_path(".grm.id");

        // "FAM_1_IND_2" should split to FID="FAM" and IID="1_IND_2"
        std::vector<std::string> ids = {"FAM_1_IND_2"};

        {
            GrmIdWriter writer(file_path);
            writer.write(ids);
        }

        auto content = read_file_content(file_path);
        REQUIRE(content == "FAM\t1_IND_2\n");
    }

    SECTION("Happy path - multiple IDs with multiple underscores")
    {
        auto file_path = files.generate_random_file_path(".grm.id");

        std::vector<std::string> ids = {
            "A_B_C_D",      // FID=A, IID=B_C_D
            "X__Y",         // FID=X, IID=_Y (double underscore)
            "TEST_1_2_3_4"  // FID=TEST, IID=1_2_3_4
        };

        {
            GrmIdWriter writer(file_path);
            writer.write(ids);
        }

        auto lines = read_file_lines(file_path);
        REQUIRE(lines.size() == 3);
        REQUIRE(lines[0] == "A\tB_C_D");
        REQUIRE(lines[1] == "X\t_Y");
        REQUIRE(lines[2] == "TEST\t1_2_3_4");
    }
}

TEST_CASE("GrmIdWriter - ID with no underscore", "[grm_id_writer][split]")
{
    FileFixture files;

    SECTION("Happy path - no underscore uses same value for FID and IID")
    {
        auto file_path = files.generate_random_file_path(".grm.id");

        std::vector<std::string> ids = {"SAMPLE123"};

        {
            GrmIdWriter writer(file_path);
            writer.write(ids);
        }

        auto content = read_file_content(file_path);
        // No underscore: both FID and IID are set to the original id
        REQUIRE(content == "SAMPLE123\tSAMPLE123\n");
    }

    SECTION("Happy path - mixed IDs with and without underscores")
    {
        auto file_path = files.generate_random_file_path(".grm.id");

        std::vector<std::string> ids = {
            "FAM1_IND1",  // Has underscore
            "NOSPLIT",    // No underscore
            "FAM2_IND2"   // Has underscore
        };

        {
            GrmIdWriter writer(file_path);
            writer.write(ids);
        }

        auto lines = read_file_lines(file_path);
        REQUIRE(lines.size() == 3);
        REQUIRE(lines[0] == "FAM1\tIND1");
        REQUIRE(lines[1] == "NOSPLIT\tNOSPLIT");
        REQUIRE(lines[2] == "FAM2\tIND2");
    }
}

TEST_CASE("GrmIdWriter - ID with leading underscore", "[grm_id_writer][split]")
{
    FileFixture files;

    SECTION("Happy path - leading underscore gives empty FID")
    {
        auto file_path = files.generate_random_file_path(".grm.id");

        std::vector<std::string> ids = {"_IND1"};

        {
            GrmIdWriter writer(file_path);
            writer.write(ids);
        }

        auto content = read_file_content(file_path);
        // FID is empty string, IID is "IND1"
        REQUIRE(content == "\tIND1\n");
    }
}

TEST_CASE("GrmIdWriter - ID with trailing underscore", "[grm_id_writer][split]")
{
    FileFixture files;

    SECTION("Happy path - trailing underscore gives empty IID")
    {
        auto file_path = files.generate_random_file_path(".grm.id");

        std::vector<std::string> ids = {"FAM1_"};

        {
            GrmIdWriter writer(file_path);
            writer.write(ids);
        }

        auto content = read_file_content(file_path);
        // FID is "FAM1", IID is empty string
        REQUIRE(content == "FAM1\t\n");
    }
}

// ============================================================================
// Edge case tests
// ============================================================================

TEST_CASE("GrmIdWriter - Empty string ID", "[grm_id_writer][edge]")
{
    FileFixture files;

    SECTION("Happy path - empty string ID")
    {
        auto file_path = files.generate_random_file_path(".grm.id");

        std::vector<std::string> ids = {""};

        {
            GrmIdWriter writer(file_path);
            writer.write(ids);
        }

        auto content = read_file_content(file_path);
        // Empty string has no underscore, so both FID and IID are empty
        REQUIRE(content == "\t\n");
    }
}

TEST_CASE("GrmIdWriter - ID with only underscore", "[grm_id_writer][edge]")
{
    FileFixture files;

    SECTION("Happy path - single underscore only")
    {
        auto file_path = files.generate_random_file_path(".grm.id");

        std::vector<std::string> ids = {"_"};

        {
            GrmIdWriter writer(file_path);
            writer.write(ids);
        }

        auto content = read_file_content(file_path);
        // Split "_" -> FID="", IID=""
        REQUIRE(content == "\t\n");
    }
}

// ============================================================================
// Output format verification tests
// ============================================================================

TEST_CASE("GrmIdWriter - Output format verification", "[grm_id_writer][format]")
{
    FileFixture files;

    SECTION("Happy path - verify tab separator and newline")
    {
        auto file_path = files.generate_random_file_path(".grm.id");

        std::vector<std::string> ids = {"A_B", "C_D"};

        {
            GrmIdWriter writer(file_path);
            writer.write(ids);
        }

        auto content = read_file_content(file_path);

        // Verify exact format: "FID\tIID\n" for each line
        REQUIRE(content == "A\tB\nC\tD\n");

        // Count tabs and newlines
        auto tab_count = std::count(content.begin(), content.end(), '\t');
        auto newline_count = std::count(content.begin(), content.end(), '\n');

        REQUIRE(tab_count == 2);
        REQUIRE(newline_count == 2);
    }

    SECTION("Happy path - verify no trailing content after last newline")
    {
        auto file_path = files.generate_random_file_path(".grm.id");

        std::vector<std::string> ids = {"X_Y"};

        {
            GrmIdWriter writer(file_path);
            writer.write(ids);
        }

        auto content = read_file_content(file_path);

        // File should end with newline
        REQUIRE(!content.empty());
        REQUIRE(content.back() == '\n');
    }
}

// ============================================================================
// Multiple write calls tests
// ============================================================================

TEST_CASE("GrmIdWriter - Multiple write calls", "[grm_id_writer][multiple]")
{
    FileFixture files;

    SECTION("Happy path - multiple write calls append content")
    {
        auto file_path = files.generate_random_file_path(".grm.id");

        std::vector<std::string> ids1 = {"FAM1_IND1", "FAM2_IND2"};
        std::vector<std::string> ids2 = {"FAM3_IND3"};
        std::vector<std::string> ids3 = {"FAM4_IND4", "FAM5_IND5"};

        {
            GrmIdWriter writer(file_path);
            writer.write(ids1);
            writer.write(ids2);
            writer.write(ids3);
        }

        auto lines = read_file_lines(file_path);
        REQUIRE(lines.size() == 5);
        REQUIRE(lines[0] == "FAM1\tIND1");
        REQUIRE(lines[1] == "FAM2\tIND2");
        REQUIRE(lines[2] == "FAM3\tIND3");
        REQUIRE(lines[3] == "FAM4\tIND4");
        REQUIRE(lines[4] == "FAM5\tIND5");
    }

    SECTION("Happy path - write empty then non-empty")
    {
        auto file_path = files.generate_random_file_path(".grm.id");

        std::vector<std::string> empty_ids;
        std::vector<std::string> ids = {"A_B"};

        {
            GrmIdWriter writer(file_path);
            writer.write(empty_ids);
            writer.write(ids);
        }

        auto content = read_file_content(file_path);
        REQUIRE(content == "A\tB\n");
    }
}

// ============================================================================
// Large input tests
// ============================================================================

TEST_CASE("GrmIdWriter - Large number of IDs", "[grm_id_writer][large]")
{
    FileFixture files;

    SECTION("Happy path - write 1000 IDs")
    {
        auto file_path = files.generate_random_file_path(".grm.id");

        constexpr size_t num_ids = 1000;
        std::vector<std::string> ids;
        ids.reserve(num_ids);

        for (size_t i = 0; i < num_ids; ++i)
        {
            ids.push_back(
                "FAM" + std::to_string(i) + "_IND" + std::to_string(i));
        }

        REQUIRE_NOTHROW(
            [&]()
            {
                GrmIdWriter writer(file_path);
                writer.write(ids);
            }());

        auto lines = read_file_lines(file_path);
        REQUIRE(lines.size() == num_ids);

        // Verify first and last lines
        REQUIRE(lines[0] == "FAM0\tIND0");
        REQUIRE(lines[num_ids - 1] == "FAM999\tIND999");
    }
}

// ============================================================================
// Special character tests
// ============================================================================

TEST_CASE(
    "GrmIdWriter - IDs with special characters",
    "[grm_id_writer][special]")
{
    FileFixture files;

    SECTION("Happy path - IDs with numbers")
    {
        auto file_path = files.generate_random_file_path(".grm.id");

        std::vector<std::string> ids = {"123_456", "FAM01_IND02"};

        {
            GrmIdWriter writer(file_path);
            writer.write(ids);
        }

        auto lines = read_file_lines(file_path);
        REQUIRE(lines.size() == 2);
        REQUIRE(lines[0] == "123\t456");
        REQUIRE(lines[1] == "FAM01\tIND02");
    }

    SECTION("Happy path - IDs with dots and dashes")
    {
        auto file_path = files.generate_random_file_path(".grm.id");

        std::vector<std::string> ids = {"FAM.1_IND-1", "A-B.C_D.E-F"};

        {
            GrmIdWriter writer(file_path);
            writer.write(ids);
        }

        auto lines = read_file_lines(file_path);
        REQUIRE(lines.size() == 2);
        REQUIRE(lines[0] == "FAM.1\tIND-1");
        REQUIRE(lines[1] == "A-B.C\tD.E-F");
    }
}
