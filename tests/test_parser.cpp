#include <string>
#include <vector>

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include "../src/data/parser.h"

using Catch::Matchers::WithinAbs;

TEST_CASE("try_parse_double function", "[parser]")
{
    SECTION("Valid double values")
    {
        REQUIRE_THAT(
            gelex::detail::try_parse_double("3.14").value(),
            WithinAbs(3.14, 1e-10));
        REQUIRE_THAT(
            gelex::detail::try_parse_double("-2.5").value(),
            WithinAbs(-2.5, 1e-10));
        REQUIRE_THAT(
            gelex::detail::try_parse_double("0").value(),
            WithinAbs(0.0, 1e-10));
        REQUIRE_THAT(
            gelex::detail::try_parse_double("1e6").value(),
            WithinAbs(1e6, 1e-10));
    }

    SECTION("Invalid double values")
    {
        REQUIRE(gelex::detail::try_parse_double("abc").has_value() == false);
        REQUIRE(
            gelex::detail::try_parse_double("3.14abc").has_value() == false);
        REQUIRE(gelex::detail::try_parse_double("").has_value() == false);
    }
}

TEST_CASE("parse_nth_double function", "[parser]")
{
    const std::string_view line = "FID\tIID\t3.14\t-2.5\t0\t1e6";

    SECTION("Valid column indices")
    {
        REQUIRE_THAT(
            gelex::detail::parse_nth_double(line, 2).value(),
            WithinAbs(3.14, 1e-10));
        REQUIRE_THAT(
            gelex::detail::parse_nth_double(line, 3).value(),
            WithinAbs(-2.5, 1e-10));
        REQUIRE_THAT(
            gelex::detail::parse_nth_double(line, 4).value(),
            WithinAbs(0.0, 1e-10));
        REQUIRE_THAT(
            gelex::detail::parse_nth_double(line, 5).value(),
            WithinAbs(1e6, 1e-10));
    }

    SECTION("Invalid column indices")
    {
        REQUIRE(gelex::detail::parse_nth_double(line, 0).has_value() == false);
        REQUIRE(gelex::detail::parse_nth_double(line, 1).has_value() == false);
        REQUIRE(gelex::detail::parse_nth_double(line, 10).has_value() == false);
    }

    SECTION("Custom delimiters")
    {
        const std::string_view comma_line = "FID,IID,3.14,-2.5,0,1e6";
        REQUIRE_THAT(
            gelex::detail::parse_nth_double(comma_line, 2, ",").value(),
            WithinAbs(3.14, 1e-10));
    }
}

TEST_CASE("parse_id function", "[parser]")
{
    const std::string_view line = "FAM001\tIND001\t3.14\t-2.5";

    SECTION("IID only mode")
    {
        auto result = gelex::detail::parse_id(line, true);
        REQUIRE(result.value() == "IND001");
    }

    SECTION("Full ID mode")
    {
        auto result = gelex::detail::parse_id(line, false);
        REQUIRE(result.value() == "FAM001_IND001");
    }

    SECTION("Insufficient columns")
    {
        const std::string_view short_line = "FAM001";
        auto result = gelex::detail::parse_id(short_line, false);
        REQUIRE(result.has_value() == false);
    }

    SECTION("Custom delimiters")
    {
        const std::string_view comma_line = "FAM001,IND001,3.14,-2.5";
        auto result = gelex::detail::parse_id(comma_line, false, ",");
        REQUIRE(result.value() == "FAM001_IND001");
    }
}

TEST_CASE("parse_header function", "[parser]")
{
    SECTION("Valid header")
    {
        const std::string_view header_line = "FID\tIID\tTrait1\tTrait2";
        auto result = gelex::detail::parse_header(header_line);

        REQUIRE(result.has_value());
        REQUIRE(result->size() == 4);
        REQUIRE(result->at(0) == "FID");
        REQUIRE(result->at(1) == "IID");
        REQUIRE(result->at(2) == "Trait1");
        REQUIRE(result->at(3) == "Trait2");
    }

    SECTION("Empty header")
    {
        const std::string_view empty_line;
        auto result = gelex::detail::parse_header(empty_line);
        REQUIRE(result.has_value() == false);
    }

    SECTION("Insufficient columns")
    {
        const std::string_view short_line = "FID";
        auto result = gelex::detail::parse_header(short_line);
        REQUIRE(result.has_value() == false);
    }

    SECTION("Wrong column names")
    {
        const std::string_view wrong_line = "ID\tNAME\tTrait1";
        auto result = gelex::detail::parse_header(wrong_line);
        REQUIRE(result.has_value() == false);
    }

    SECTION("Custom delimiters")
    {
        const std::string_view comma_line = "FID,IID,Trait1,Trait2";
        auto result = gelex::detail::parse_header(comma_line, ",");
        REQUIRE(result.has_value());
        REQUIRE(result->size() == 4);
    }
}

TEST_CASE("count_num_columns function", "[parser]")
{
    SECTION("Basic counting")
    {
        REQUIRE(gelex::detail::count_num_columns("FID\tIID\tTrait1") == 3);
        REQUIRE(gelex::detail::count_num_columns("a\tb\tc\td\te") == 5);
        REQUIRE(gelex::detail::count_num_columns("single") == 1);
    }

    SECTION("Empty and whitespace handling")
    {
        REQUIRE(gelex::detail::count_num_columns("") == 0);
        REQUIRE(gelex::detail::count_num_columns("\t\t") == 0);
        REQUIRE(gelex::detail::count_num_columns("a\t\tb") == 2);
    }

    SECTION("Custom delimiters")
    {
        REQUIRE(gelex::detail::count_num_columns("a,b,c", ",") == 3);
        REQUIRE(gelex::detail::count_num_columns("a;b;c", ";") == 3);
    }
}

TEST_CASE("parse_string function", "[parser]")
{
    const std::string_view line = "FID\tIID\t3.14\t-2.5\t0";

    SECTION("Full line parsing")
    {
        auto result = gelex::detail::parse_string(line);
        REQUIRE(result.size() == 5);
        REQUIRE(result[0] == "FID");
        REQUIRE(result[1] == "IID");
        REQUIRE(result[2] == "3.14");
        REQUIRE(result[3] == "-2.5");
        REQUIRE(result[4] == "0");
    }

    SECTION("With column offset")
    {
        auto result = gelex::detail::parse_string(line, 2);
        REQUIRE(result.size() == 3);
        REQUIRE(result[0] == "3.14");
        REQUIRE(result[1] == "-2.5");
        REQUIRE(result[2] == "0");
    }

    SECTION("Custom delimiters")
    {
        const std::string_view comma_line = "FID,IID,3.14,-2.5,0";
        auto result = gelex::detail::parse_string(comma_line, 0, ",");
        REQUIRE(result.size() == 5);
        REQUIRE(result[0] == "FID");
    }

    SECTION("Empty and whitespace handling")
    {
        auto empty_result = gelex::detail::parse_string("");
        REQUIRE(empty_result.empty());

        auto whitespace_result = gelex::detail::parse_string("\t\t");
        REQUIRE(whitespace_result.empty());
    }
}

TEST_CASE("parse_all_doubles function", "[parser]")
{
    const std::string_view line = "FID\tIID\t3.14\t-2.5\t0\t1e6";

    SECTION("Valid double parsing")
    {
        auto result = gelex::detail::parse_all_doubles(line, 2);
        REQUIRE(result.has_value());
        REQUIRE(result->size() == 4);
        REQUIRE_THAT(result->at(0), WithinAbs(3.14, 1e-10));
        REQUIRE_THAT(result->at(1), WithinAbs(-2.5, 1e-10));
        REQUIRE_THAT(result->at(2), WithinAbs(0.0, 1e-10));
        REQUIRE_THAT(result->at(3), WithinAbs(1e6, 1e-10));
    }

    SECTION("Invalid double values")
    {
        const std::string_view bad_line = "FID\tIID\t3.14\tabc\t0";
        auto result = gelex::detail::parse_all_doubles(bad_line, 2);
        REQUIRE(result.has_value() == false);
    }

    SECTION("Empty input")
    {
        auto result = gelex::detail::parse_all_doubles("", 0);
        REQUIRE(result.has_value());
        REQUIRE(result->empty());
    }

    SECTION("Custom delimiters")
    {
        const std::string_view comma_line = "FID,IID,3.14,-2.5,0,1e6";
        auto result = gelex::detail::parse_all_doubles(comma_line, 2, ",");
        REQUIRE(result.has_value());
        REQUIRE(result->size() == 4);
    }

    SECTION("Error message includes column index")
    {
        const std::string_view bad_line = "FID\tIID\t3.14\tabc\t0";
        auto result = gelex::detail::parse_all_doubles(bad_line, 2);
        REQUIRE(result.has_value() == false);
        REQUIRE(result.error().message.find("at column") != std::string::npos);
    }
}
