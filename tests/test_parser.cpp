#include <cmath>
#include <filesystem>
#include <fstream>
#include <string_view>
#include <vector>

#include <catch2/catch_test_macros.hpp>

#include "../src/data/parser.h"
#include "file_fixture.h"
#include "gelex/exception.h"

namespace fs = std::filesystem;

using namespace gelex::detail;  // NOLINT
using gelex::test::FileFixture;

TEST_CASE("File Stream I/O", "[parser]")
{
    FileFixture files;
    SECTION("Happy path - open existing file")
    {
        auto file_path = files.create_text_file("content");

        REQUIRE_NOTHROW(
            [&]()
            {
                auto file = open_file<std::ifstream>(file_path, std::ios::in);
                REQUIRE(file.is_open());
            }());
    }

    SECTION("Happy path - open file with custom buffer")
    {
        auto file_path = files.create_text_file("content");
        std::vector<char> buffer(4096);  // 4KB buffer

        REQUIRE_NOTHROW(
            [&]()
            {
                auto file
                    = open_file<std::ifstream>(file_path, std::ios::in, buffer);
                REQUIRE(file.is_open());
            }());
    }

    SECTION("Happy path - open file for writing")
    {
        auto file_path = files.get_test_dir() / "output.txt";

        REQUIRE_NOTHROW(
            [&]()
            {
                auto file = open_file<std::ofstream>(file_path, std::ios::out);
                REQUIRE(file.is_open());
                file << "test content";
            }());

        REQUIRE(std::filesystem::exists(file_path));
        REQUIRE(std::filesystem::file_size(file_path) > 0);
    }

    SECTION("Happy path - open file for appending")
    {
        auto file_path = files.create_text_file("initial content");

        REQUIRE_NOTHROW(
            [&]()
            {
                auto file = open_file<std::ofstream>(
                    file_path, std::ios::out | std::ios::app);
                REQUIRE(file.is_open());
                file << "\nappended content";
            }());

        // Read back to verify
        std::ifstream check(file_path);
        std::string content(
            (std::istreambuf_iterator<char>(check)),
            std::istreambuf_iterator<char>());
        REQUIRE(content.find("appended content") != std::string::npos);
    }

    SECTION("Happy path - open file for reading and writing")
    {
        auto file_path = files.create_text_file("initial content");

        REQUIRE_NOTHROW(
            [&]()
            {
                auto file = open_file<std::fstream>(
                    file_path, std::ios::in | std::ios::out);
                REQUIRE(file.is_open());

                // Read existing content
                std::string line;
                std::getline(file, line);
                REQUIRE(line == "initial content");

                // Write new content
                file.seekp(0, std::ios::end);
                file << "\nnew content";
            }());
    }

    SECTION("Exception - file not found")
    {
        REQUIRE_THROWS_AS(
            open_file<std::ifstream>("non_existent_path", std::ios::in),
            gelex::FileNotFoundException);
    }

    SECTION("Exception - empty file")
    {
        auto empty_file_path = files.create_empty_file();

        REQUIRE_THROWS_AS(
            open_file<std::ifstream>(empty_file_path, std::ios::in),
            gelex::FileFormatException);
    }

    SECTION("Exception - directory instead of file")
    {
        auto dir_path = files.get_test_dir();

        REQUIRE_THROWS_AS(
            open_file<std::ifstream>(dir_path, std::ios::in),
            gelex::FileOpenException);
        REQUIRE_THROWS_AS(
            open_file<std::ofstream>(dir_path, std::ios::out),
            gelex::FileOpenException);
    }

    SECTION("Edge case - empty buffer")
    {
        auto file_path = files.create_text_file("content");
        std::vector<char> empty_buffer;

        REQUIRE_NOTHROW(
            [&]()
            {
                auto file = open_file<std::ifstream>(
                    file_path, std::ios::in, empty_buffer);
                REQUIRE(file.is_open());
            }());
    }

    SECTION("Edge case - writing to empty file is allowed")
    {
        auto empty_file_path = files.create_empty_file();

        REQUIRE_NOTHROW(
            [&]()
            {
                auto file
                    = open_file<std::ofstream>(empty_file_path, std::ios::out);
                REQUIRE(file.is_open());
                file << "writing to empty file";
            }());
    }
}

TEST_CASE("Parser Line Counting Tests", "[parser]")
{
    FileFixture files;
    SECTION("Happy path - count lines")
    {
        auto file_path = files.create_text_file("line1\nline2\nline3\nline4\n");
        REQUIRE(count_total_lines(file_path) == 4);
    }

    SECTION("Happy path - file without trailing newline")
    {
        auto file_path = files.create_text_file("line1\nline2\nline3");
        REQUIRE(count_total_lines(file_path) == 3);
    }
}

TEST_CASE("Parser Column Count Tests", "[parser]")
{
    SECTION("Happy path - tab delimiter")
    {
        std::string_view line = "FID\tIID\tPhenotype\tCovariate";

        REQUIRE(count_num_columns(line, '\t') == 4);
    }

    SECTION("Happy path - comma delimiter")
    {
        std::string_view line = "FID,IID,Phenotype,Covariate";

        REQUIRE(count_num_columns(line, ',') == 4);
    }

    SECTION("Edge case - empty line")
    {
        std::string_view line;

        REQUIRE(count_num_columns(line, '\t') == 0);
    }

    SECTION("Edge case - trailing delimiter")
    {
        std::string_view line = "FID\tIID\tPhenotype\t";

        REQUIRE(count_num_columns(line, '\t') == 4);
    }
}

TEST_CASE("Parser Double Parsing Tests", "[parser]")
{
    SECTION("Happy path - valid numbers")
    {
        REQUIRE(parse_number("42") == 42.0);
        REQUIRE(parse_number("3.14") == 3.14);
        REQUIRE(parse_number("-1.5") == -1.5);
        REQUIRE(parse_number("1.23e-4") == 1.23e-4);
        REQUIRE(std::isnan(parse_number("nan")));
    }

    SECTION("Exception - invalid numbers")
    {
        REQUIRE_THROWS_AS(parse_number("abc"), gelex::NumberParseException);
        REQUIRE_THROWS_AS(parse_number("1.2.3"), gelex::NumberParseException);
        REQUIRE_THROWS_AS(parse_number(""), gelex::NumberParseException);
        REQUIRE_THROWS_AS(parse_number(" "), gelex::NumberParseException);
    }
}

TEST_CASE("Parser Column Double Parsing Tests", "[parser]")
{
    SECTION("Happy path - extract valid doubles")
    {
        std::string_view line = "FID\tIID\t2.5\t1.0\t0.5";

        REQUIRE(parse_nth_double(line, 2, '\t') == 2.5);
        REQUIRE(parse_nth_double(line, 3, '\t') == 1.0);
        REQUIRE(parse_nth_double(line, 4, '\t') == 0.5);
    }

    SECTION("Happy path - custom delimiter")
    {
        std::string_view line = "FID,IID,2.5,1.0,0.5";

        REQUIRE(parse_nth_double(line, 2, ',') == 2.5);
    }

    SECTION("Exception - column index out of range")
    {
        std::string_view line = "FID\tIID\t2.5";

        REQUIRE_THROWS_AS(
            parse_nth_double(line, 5, '\t'), gelex::ColumnRangeException);
    }

    SECTION("Exception - invalid number at column")
    {
        std::string_view line = "FID\tIID\tinvalid";

        REQUIRE_THROWS_AS(
            parse_nth_double(line, 2, '\t'), gelex::NumberParseException);
    }
}

TEST_CASE("Parser All Doubles Parsing Tests", "[parser]")
{
    SECTION("Happy path - parse all doubles")
    {
        std::string_view line = "2.5\t1.0\t0.5\t-1.2";
        std::vector<double> doubles;
        parse_all_doubles(line, doubles, 0, '\t');
        REQUIRE(doubles.size() == 4);
        REQUIRE(doubles[0] == 2.5);
        REQUIRE(doubles[1] == 1.0);
        REQUIRE(doubles[2] == 0.5);
        REQUIRE(doubles[3] == -1.2);
    }

    SECTION("Happy path - parse with column offset")
    {
        std::string_view line = "FID\tIID\t2.5\t1.0\t0.5";
        std::vector<double> doubles;
        parse_all_doubles(line, doubles, 2, '\t');
        REQUIRE(doubles.size() == 3);
        REQUIRE(doubles[0] == 2.5);
        REQUIRE(doubles[1] == 1.0);
        REQUIRE(doubles[2] == 0.5);
    }

    SECTION("Exception - invalid number in any column")
    {
        std::string_view line = "2.5\tinvalid\t0.5";

        std::vector<double> doubles;
        REQUIRE_THROWS_AS(
            parse_all_doubles(line, doubles, 0, '\t'),
            gelex::DataParseException);
    }

    SECTION("Edge case - empty line")
    {
        std::string_view line;

        std::vector<double> doubles;
        parse_all_doubles(line, doubles, 0, '\t');
        REQUIRE(doubles.empty());
    }
}

TEST_CASE("Parser Header Tests", "[parser]")
{
    SECTION("Happy path - valid header")
    {
        std::string_view line = "FID\tIID\tPhenotype\tCovariate";

        auto header = parse_header(line, '\t');
        REQUIRE(header.size() == 4);
        REQUIRE(header[0] == "FID");
        REQUIRE(header[1] == "IID");
        REQUIRE(header[2] == "Phenotype");
        REQUIRE(header[3] == "Covariate");
    }

    SECTION("Happy path - custom delimiter")
    {
        std::string_view line = "FID,IID,Phenotype,Covariate";

        auto header = parse_header(line, ',');
        REQUIRE(header.size() == 4);
        REQUIRE(header[0] == "FID");
        REQUIRE(header[1] == "IID");
    }

    SECTION("Exception - missing FID column")
    {
        std::string_view line = "ID\tIID\tPhenotype";

        REQUIRE_THROWS_AS(
            parse_header(line, '\t'), gelex::HeaderFormatException);
    }

    SECTION("Exception - missing IID column")
    {
        std::string_view line = "FID\tSample\tPhenotype";

        REQUIRE_THROWS_AS(
            parse_header(line, '\t'), gelex::HeaderFormatException);
    }

    SECTION("Exception - insufficient columns")
    {
        std::string_view line = "FID";

        REQUIRE_THROWS_AS(
            parse_header(line, '\t'), gelex::HeaderFormatException);
    }

    SECTION("Exception - empty header")
    {
        std::string_view line;

        REQUIRE_THROWS_AS(
            parse_header(line, '\t'), gelex::HeaderFormatException);
    }
}

TEST_CASE("Parser ID Parsing Tests", "[parser]")
{
    SECTION("Happy path - full ID with tab delimiter")
    {
        std::string_view line = "1\t2\t2.5\t1.0";

        REQUIRE(parse_id(line, false, '\t') == "1_2");
    }

    SECTION("Happy path - IID only with tab delimiter")
    {
        std::string_view line = "1\t2\t2.5\t1.0";

        REQUIRE(parse_id(line, true, '\t') == "2");
    }

    SECTION("Happy path - custom delimiter")
    {
        std::string_view line = "1,2,2.5,1.0";

        REQUIRE(parse_id(line, false, ',') == "1_2");
    }

    SECTION("Exception - insufficient columns")
    {
        std::string_view line = "1";

        REQUIRE_THROWS_AS(
            parse_id(line, false, '\t'), gelex::FileFormatException);
    }
}

TEST_CASE("Parser String Parsing Tests", "[parser]")
{
    SECTION("Happy path - parse all columns")
    {
        std::string_view line = "FID\tIID\tPhenotype\tCovariate";
        std::vector<std::string_view> strings;
        parse_string(line, strings, 0, '\t');

        REQUIRE(strings.size() == 4);
        REQUIRE(strings[0] == "FID");
        REQUIRE(strings[1] == "IID");
        REQUIRE(strings[2] == "Phenotype");
        REQUIRE(strings[3] == "Covariate");
    }

    SECTION("Happy path - parse with column offset")
    {
        std::string_view line = "FID\tIID\tPhenotype\tCovariate";
        std::vector<std::string_view> strings;
        parse_string(line, strings, 2, '\t');

        REQUIRE(strings.size() == 2);
        REQUIRE(strings[0] == "Phenotype");
        REQUIRE(strings[1] == "Covariate");
    }

    SECTION("Edge case - offset beyond available columns")
    {
        std::string_view line = "FID\tIID\tPhenotype";
        std::vector<std::string_view> strings;
        parse_string(line, strings, 5, '\t');

        REQUIRE(strings.empty());
    }

    SECTION("Edge case - empty columns")
    {
        std::string_view line = "FID\t\tPhenotype\t";
        std::vector<std::string_view> strings;

        REQUIRE_THROWS_AS(
            parse_string(line, strings, 0, '\t'), gelex::DataParseException);
    }
}
