#include <string_view>
#include <vector>

#include <catch2/catch_test_macros.hpp>

#include "../src/data/parser.h"

TEST_CASE("try_parse_double tests", "[parser]")
{
    SECTION("Happy Path: 输入是合法的双精度浮点数字符串")
    {
        // 目的：验证函数能否正确解析包含有效浮点数的字符串。
        // 预期：返回一个包含解析后 double 值的 std::expected 对象。
        std::string_view valid_double_str = "123.45";
        auto result = gelex::detail::try_parse_double(valid_double_str);
        REQUIRE(result.has_value());
        CHECK(result.value() == 123.45);
    }

    SECTION("错误场景一: 输入是纯粹的非数字字符串")
    {
        // 目的：验证函数如何处理完全不包含数字的字符串。
        // 预期：返回一个包含 ParseError::NotNumber 错误的 std::expected 对象。
        std::string_view invalid_str = "not_a_number";
        auto result = gelex::detail::try_parse_double(invalid_str);
        REQUIRE_FALSE(result.has_value());
        CHECK(result.error() == gelex::detail::ParseError::NotNumber);
    }

    SECTION("错误场景二: 输入的字符串包含数字和尾随的非数字字符")
    {
        // 目的：验证函数是否能拒绝部分合法的输入，确保整个字符串被完整解析。
        // 预期：返回一个包含 ParseError::NotNumber 错误的 std::expected
        // 对象，因为整个 token 并非一个有效的数字。
        std::string_view mixed_str = "123.45abc";
        auto result = gelex::detail::try_parse_double(mixed_str);
        REQUIRE_FALSE(result.has_value());
        CHECK(result.error() == gelex::detail::ParseError::NotNumber);
    }
}

TEST_CASE("parse_nth_double tests", "[parser]")
{
    SECTION("Happy Path: 从有效列索引成功解析 double")
    {
        // 目的：验证函数能否从一个由分隔符分割的字符串中，根据给定的索引正确解析出浮点数。
        // 预期：返回一个包含目标索引处 double 值的 std::expected 对象。
        std::string_view line = "1.1\t2.2\t3.3";
        size_t column_index = 1;
        auto result = gelex::detail::parse_nth_double(line, column_index, "\t");
        REQUIRE(result.has_value());
        CHECK(result.value() == 2.2);
    }

    SECTION("错误场景一: 指定索引处的 token 不是一个有效的数字")
    {
        // 目的：验证当目标索引位置的子字符串无法被解析为浮点数时，函数的错误处理能力。
        // 预期：返回一个包含 ParseError::NotNumber 错误的 std::expected 对象。
        std::string_view line = "1.1\ttext\t3.3";
        size_t column_index = 1;
        auto result = gelex::detail::parse_nth_double(line, column_index, "\t");
        REQUIRE_FALSE(result.has_value());
        CHECK(result.error() == gelex::detail::ParseError::NotNumber);
    }

    SECTION("错误场景二: 列索引超出范围")
    {
        // 目的：验证当提供的索引值大于等于字符串中 token 的数量时，函数的行为。
        // 预期：返回一个包含 ParseError::InvalidColumn 错误的 std::expected
        // 对象。
        std::string_view line = "1.1\t2.2\t3.3";
        size_t column_index = 5;
        auto result = gelex::detail::parse_nth_double(line, column_index, "\t");
        REQUIRE_FALSE(result.has_value());
        CHECK(result.error() == gelex::detail::ParseError::InvalidColumn);
    }
}

TEST_CASE("parse_string tests", "[parser]")
{
    SECTION("Happy Path: 使用默认偏移量解析一行中的所有 token")
    {
        // 目的：验证函数能够正确地使用分隔符分割字符串，并返回所有非空子串的
        // vector。 预期：返回一个包含所有分割后的 std::string_view 的 vector。
        std::string_view line = "one two  three";
        std::vector<std::string_view> expected = {"one", "two", "three"};
        auto result = gelex::detail::parse_string(line, 0, " ");
        REQUIRE(result == expected);
    }

    SECTION("错误场景一: 输入行为空字符串")
    {
        // 目的：验证函数在处理空输入时的行为。由于该函数无法返回错误，此为一个边界条件测试。
        // 预期：返回一个空的 vector。
        std::string_view line = "";
        auto result = gelex::detail::parse_string(line, 0, " ");
        REQUIRE(result.empty());
    }

    SECTION("错误场景二: 列偏移量大于 token 数量")
    {
        // 目的：验证当指定的偏移量超出实际 token
        // 数量时函数的行为。此为另一个边界条件测试。 预期：返回一个空的
        // vector，因为 drop 操作会消耗所有 token。
        std::string_view line = "one two three";
        auto result = gelex::detail::parse_string(line, 5, " ");
        REQUIRE(result.empty());
    }
}

TEST_CASE("parse_all_doubles tests", "[parser]")
{
    SECTION("Happy Path: 所有 token 都是有效的 double")
    {
        // 目的：验证函数能够成功解析一个所有子串都是有效浮点数的字符串。
        // 预期：返回一个 std::expected 对象，其内部包含一个持有所有解析后
        // double 值的 vector。
        std::string_view line = "1.1\t2.2\t3.3";
        std::vector<double> expected = {1.1, 2.2, 3.3};
        auto result = gelex::detail::parse_all_doubles(line, 0, "\t");
        REQUIRE(result.has_value());
        CHECK(result.value() == expected);
    }

    SECTION("错误场景一: 行中包含非数字 token")
    {
        // 目的：验证当字符串中混有非数字子串时，函数能够立即中止并报告错误。
        // 预期：返回一个包含 ParseError::NotNumber 错误的 std::expected 对象。
        std::string_view line = "1.1\tnot-a-number\t3.3";
        auto result = gelex::detail::parse_all_doubles(line, 0, "\t");
        REQUIRE_FALSE(result.has_value());
        CHECK(result.error() == gelex::detail::ParseError::NotNumber);
    }

    SECTION("错误场景二: 行中包含格式错误的数字 token")
    {
        // 目的：验证函数能否识别并拒绝格式不正确的数字（例如，多个小数点）。
        // 预期：返回一个包含 ParseError::NotNumber 错误的 std::expected 对象。
        std::string_view line = "1.1\t1.2.3\t3.3";
        auto result = gelex::detail::parse_all_doubles(line, 0, "\t");
        REQUIRE_FALSE(result.has_value());
        CHECK(result.error() == gelex::detail::ParseError::NotNumber);
    }
}
