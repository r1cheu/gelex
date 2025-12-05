#include <vector>

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_exception.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>

#include "file_fixture.h"
#include "gelex/exception.h"

#include "../src/data/loader/snp_effect_loader.h"

namespace fs = std::filesystem;

using namespace gelex;          // For SnpEffects
using namespace gelex::detail;  // NOLINT
using Catch::Matchers::EndsWith;
using gelex::test::FileFixture;

// Helper function to create test SNP effect files
std::string create_snp_effect_content(
    const std::string& header,
    const std::vector<std::string>& rows)
{
    std::string content = header + "\n";
    for (const auto& row : rows)
    {
        content += row + "\n";
    }
    return content;
}

TEST_CASE("ColumnIndices - has_required_columns", "[data][snp_effect]")
{
    SECTION("Happy path - all required columns present")
    {
        ColumnIndices indices{
            .id = 0, .a1 = 1, .a2 = 2, .a1frq = 3, .add = 4, .dom = 5};

        REQUIRE(indices.has_required_columns() == true);
    }

    SECTION("Happy path - all required columns present, dom optional")
    {
        ColumnIndices indices{
            .id = 0, .a1 = 1, .a2 = 2, .a1frq = 3, .add = 4, .dom = -1};

        REQUIRE(indices.has_required_columns() == true);
    }

    SECTION("Exception path - missing ID column")
    {
        ColumnIndices indices{
            .id = -1, .a1 = 0, .a2 = 1, .a1frq = 2, .add = 3, .dom = 4};

        REQUIRE(indices.has_required_columns() == false);
    }

    SECTION("Exception path - missing A1 column")
    {
        ColumnIndices indices{
            .id = 0, .a1 = -1, .a2 = 1, .a1frq = 2, .add = 3, .dom = 4};

        REQUIRE(indices.has_required_columns() == false);
    }

    SECTION("Exception path - missing A2 column")
    {
        ColumnIndices indices{
            .id = 0, .a1 = 1, .a2 = -1, .a1frq = 2, .add = 3, .dom = 4};

        REQUIRE(indices.has_required_columns() == false);
    }

    SECTION("Exception path - missing A1Frq column")
    {
        ColumnIndices indices{
            .id = 0, .a1 = 1, .a2 = 2, .a1frq = -1, .add = 3, .dom = 4};

        REQUIRE(indices.has_required_columns() == false);
    }

    SECTION("Exception path - missing Add column")
    {
        ColumnIndices indices{
            .id = 0, .a1 = 1, .a2 = 2, .a1frq = 3, .add = -1, .dom = 4};

        REQUIRE(indices.has_required_columns() == false);
    }
}

TEST_CASE("ColumnIndices - max_required_index", "[data][snp_effect]")
{
    SECTION("Happy path - all columns present")
    {
        ColumnIndices indices{
            .id = 0, .a1 = 1, .a2 = 2, .a1frq = 3, .add = 4, .dom = 5};

        REQUIRE(indices.max_required_index() == 5);
    }

    SECTION("Happy path - dom column has highest index")
    {
        ColumnIndices indices{
            .id = 2, .a1 = 0, .a2 = 1, .a1frq = 3, .add = 4, .dom = 5};

        REQUIRE(indices.max_required_index() == 5);
    }

    SECTION("Happy path - add column has highest index, no dom")
    {
        ColumnIndices indices{
            .id = 0, .a1 = 1, .a2 = 2, .a1frq = 3, .add = 4, .dom = -1};

        REQUIRE(indices.max_required_index() == 4);
    }

    SECTION("Happy path - a1frq column has highest index")
    {
        ColumnIndices indices{
            .id = 0, .a1 = 1, .a2 = 2, .a1frq = 5, .add = 3, .dom = 4};

        REQUIRE(indices.max_required_index() == 5);
    }
}

TEST_CASE(
    "SnpEffectLoader - Constructor and basic loading",
    "[data][snp_effect]")
{
    FileFixture files;

    SECTION("Happy path - load complete file with all columns")
    {
        std::string content = create_snp_effect_content(
            "ID\tA1\tA2\tA1Frq\tAdd\tDom",
            {"rs001\tA\tC\t0.25\t0.123\t0.045",
             "rs002\tT\tG\t0.75\t-0.456\t0.089",
             "rs003\tC\tA\t0.50\t0.789\t-0.012"});

        auto file_path = files.create_text_file(content, ".snp.eff");

        REQUIRE_NOTHROW(
            [&]()
            {
                SnpEffectLoader loader(file_path);
                const auto& effects = loader.effects();

                REQUIRE(effects.size() == 3);
                REQUIRE(loader.has_dom_effects() == true);

                // Check first SNP
                auto it1 = effects.find("rs001");
                REQUIRE(it1 != effects.end());
                REQUIRE(it1->second.index == 0);
                REQUIRE(it1->second.A1 == 'A');
                REQUIRE(it1->second.A2 == 'C');
                REQUIRE(it1->second.A1freq == 0.25);
                REQUIRE(it1->second.add == 0.123);
                REQUIRE(it1->second.dom == 0.045);

                // Check second SNP
                auto it2 = effects.find("rs002");
                REQUIRE(it2 != effects.end());
                REQUIRE(it2->second.index == 1);
                REQUIRE(it2->second.A1 == 'T');
                REQUIRE(it2->second.A2 == 'G');
                REQUIRE(it2->second.A1freq == 0.75);
                REQUIRE(it2->second.add == -0.456);
                REQUIRE(it2->second.dom == 0.089);

                // Check third SNP
                auto it3 = effects.find("rs003");
                REQUIRE(it3 != effects.end());
                REQUIRE(it3->second.index == 2);
                REQUIRE(it3->second.A1 == 'C');
                REQUIRE(it3->second.A2 == 'A');
                REQUIRE(it3->second.A1freq == 0.50);
                REQUIRE(it3->second.add == 0.789);
                REQUIRE(it3->second.dom == -0.012);
            }());
    }

    SECTION("Happy path - load file without Dom column")
    {
        std::string content = create_snp_effect_content(
            "ID\tA1\tA2\tA1Frq\tAdd",
            {"rs101\tG\tT\t0.33\t0.111",
             "rs102\tA\tC\t0.67\t-0.222",
             "rs103\tT\tA\t0.90\t0.333"});

        auto file_path = files.create_text_file(content, ".snp.eff");

        REQUIRE_NOTHROW(
            [&]()
            {
                SnpEffectLoader loader(file_path);
                const auto& effects = loader.effects();

                REQUIRE(effects.size() == 3);
                REQUIRE(loader.has_dom_effects() == false);

                // Check first SNP
                auto it1 = effects.find("rs101");
                REQUIRE(it1 != effects.end());
                REQUIRE(it1->second.index == 0);
                REQUIRE(it1->second.A1 == 'G');
                REQUIRE(it1->second.A2 == 'T');
                REQUIRE(it1->second.A1freq == 0.33);
                REQUIRE(it1->second.add == 0.111);
                REQUIRE(std::isnan(it1->second.dom));

                // Check second SNP
                auto it2 = effects.find("rs102");
                REQUIRE(it2 != effects.end());
                REQUIRE(it2->second.index == 1);
                REQUIRE(it2->second.A1 == 'A');
                REQUIRE(it2->second.A2 == 'C');
                REQUIRE(it2->second.A1freq == 0.67);
                REQUIRE(it2->second.add == -0.222);
                REQUIRE(std::isnan(it2->second.dom));
            }());
    }

    SECTION("Happy path - take_effects moves data")
    {
        std::string content = create_snp_effect_content(
            "ID\tA1\tA2\tA1Frq\tAdd\tDom",
            {"rs301\tA\tC\t0.25\t0.123\t0.045",
             "rs302\tT\tG\t0.75\t-0.456\t0.089"});

        auto file_path = files.create_text_file(content, ".snp.eff");

        SnpEffectLoader loader(file_path);
        SnpEffects effects = std::move(loader).take_effects();

        REQUIRE(effects.size() == 2);
        REQUIRE(loader.effects().empty());  // Effects should be moved out
    }
}

TEST_CASE("SnpEffectLoader - Error handling", "[data][snp_effect]")
{
    FileFixture files;

    SECTION("Exception path - missing required column in header")
    {
        std::string content = create_snp_effect_content(
            "ID\tA1\tA2\tA1Frq",  // Missing Add column
            {"rs401\tA\tC\t0.25", "rs402\tT\tG\t0.75"});

        auto file_path = files.create_text_file(content, ".snp.eff");

        REQUIRE_THROWS_MATCHES(
            [&]() { SnpEffectLoader loader(file_path); }(),
            gelex::FileFormatException,
            Catch::Matchers::MessageMatches(
                EndsWith("missing required columns (ID, A1, A2, A1Frq, Add)")));
    }

    SECTION("Exception path - insufficient columns in data row")
    {
        std::string content = create_snp_effect_content(
            "ID\tA1\tA2\tA1Frq\tAdd\tDom",
            {"rs501\tA\tC\t0.25\t0.123\t0.045",
             "rs502\tT\tG\t0.75\t-0.456",  // Missing Dom column
             "rs503\tC\tA\t0.50\t0.789\t-0.012"});

        auto file_path = files.create_text_file(content, ".snp.eff");

        REQUIRE_THROWS_MATCHES(
            [&]() { SnpEffectLoader loader(file_path); }(),
            gelex::FileFormatException,
            Catch::Matchers::MessageMatches(EndsWith(
                "has insufficient columns. Expected at least 6, got 5")));
    }

    SECTION("Exception path - invalid number format in A1Frq")
    {
        std::string content = create_snp_effect_content(
            "ID\tA1\tA2\tA1Frq\tAdd\tDom",
            {"rs601\tA\tC\t0.25\t0.123\t0.045",
             "rs602\tT\tG\tinvalid\t-0.456\t0.089",  // Invalid A1Frq
             "rs603\tC\tA\t0.50\t0.789\t-0.012"});

        auto file_path = files.create_text_file(content, ".snp.eff");

        REQUIRE_THROWS_MATCHES(
            [&]() { SnpEffectLoader loader(file_path); }(),
            gelex::FileFormatException,
            Catch::Matchers::MessageMatches(EndsWith("as number")));
    }

    SECTION("Exception path - invalid number format in Add")
    {
        std::string content = create_snp_effect_content(
            "ID\tA1\tA2\tA1Frq\tAdd\tDom",
            {"rs701\tA\tC\t0.25\t0.123\t0.045",
             "rs702\tT\tG\t0.75\tnot_a_number\t0.089",  // Invalid Add
             "rs703\tC\tA\t0.50\t0.789\t-0.012"});

        auto file_path = files.create_text_file(content, ".snp.eff");

        REQUIRE_THROWS_MATCHES(
            [&]() { SnpEffectLoader loader(file_path); }(),
            gelex::FileFormatException,
            Catch::Matchers::MessageMatches(EndsWith("as number")));
    }

    SECTION("Exception path - invalid number format in Dom")
    {
        std::string content = create_snp_effect_content(
            "ID\tA1\tA2\tA1Frq\tAdd\tDom",
            {"rs801\tA\tC\t0.25\t0.123\t0.045",
             "rs802\tT\tG\t0.75\t-0.456\tinvalid",  // Invalid Dom
             "rs803\tC\tA\t0.50\t0.789\t-0.012"});

        auto file_path = files.create_text_file(content, ".snp.eff");

        REQUIRE_THROWS_MATCHES(
            [&]() { SnpEffectLoader loader(file_path); }(),
            gelex::FileFormatException,
            Catch::Matchers::MessageMatches(EndsWith("as number")));
    }

    SECTION("Exception path - empty file")
    {
        auto file_path = files.create_empty_file(".snp.eff");

        REQUIRE_THROWS_MATCHES(
            [&]() { SnpEffectLoader loader(file_path); }(),
            gelex::FileFormatException,
            Catch::Matchers::MessageMatches(EndsWith("is empty")));
    }

    SECTION("Exception path - file with only header")
    {
        std::string content = "ID\tA1\tA2\tA1Frq\tAdd\tDom\n";
        auto file_path = files.create_text_file(content, ".snp.eff");

        REQUIRE_NOTHROW(
            [&]()
            {
                SnpEffectLoader loader(file_path);
                REQUIRE(loader.effects().empty());
            }());
    }
}

TEST_CASE("SnpEffectLoader - Column order variations", "[data][snp_effect]")
{
    FileFixture files;

    SECTION("Happy path - different column order")
    {
        std::string content = create_snp_effect_content(
            "A1Frq\tAdd\tID\tA2\tA1\tDom",  // Different order
            {"0.25\t0.123\trs1001\tC\tA\t0.045",
             "0.75\t-0.456\trs1002\tG\tT\t0.089",
             "0.50\t0.789\trs1003\tA\tC\t-0.012"});

        auto file_path = files.create_text_file(content, ".snp.eff");

        REQUIRE_NOTHROW(
            [&]()
            {
                SnpEffectLoader loader(file_path);
                const auto& effects = loader.effects();

                REQUIRE(effects.size() == 3);
                REQUIRE(loader.has_dom_effects() == true);

                // Check first SNP
                auto it1 = effects.find("rs1001");
                REQUIRE(it1 != effects.end());
                REQUIRE(it1->second.A1 == 'A');
                REQUIRE(it1->second.A2 == 'C');
                REQUIRE(it1->second.A1freq == 0.25);
                REQUIRE(it1->second.add == 0.123);
                REQUIRE(it1->second.dom == 0.045);
            }());
    }

    SECTION("Happy path - extra columns ignored")
    {
        std::string content = create_snp_effect_content(
            "ID\tA1\tA2\tA1Frq\tAdd\tExtra3\tDom\tExtra1\tExtra2",  // Extra
                                                                    // columns
            {"rs1201\tA\tC\t0.25\t0.123\t0.03\t0.045\tignore1\tignore2",
             "rs1202\tT\tG\t0.75\t-0.456\t0.02\t0.089\tignore3\tignore4",
             "rs1203\tC\tA\t0.50\t0.789\t0.03\t-0.012\tignore5\tignore6"});

        auto file_path = files.create_text_file(content, ".snp.eff");

        REQUIRE_NOTHROW(
            [&]()
            {
                SnpEffectLoader loader(file_path);
                const auto& effects = loader.effects();

                REQUIRE(effects.size() == 3);

                auto it1 = effects.find("rs1201");
                REQUIRE(it1 != effects.end());
                REQUIRE(it1->second.A1 == 'A');
                REQUIRE(it1->second.A2 == 'C');
                REQUIRE(it1->second.A1freq == 0.25);
                REQUIRE(it1->second.add == 0.123);
                REQUIRE(it1->second.dom == 0.045);
            }());
    }
}

TEST_CASE("SnpEffectLoader - Edge cases", "[data][snp_effect]")
{
    FileFixture files;

    SECTION("Happy path - single SNP")
    {
        std::string content = create_snp_effect_content(
            "ID\tA1\tA2\tA1Frq\tAdd\tDom",
            {"rs1301\tA\tC\t0.25\t0.123\t0.045"});

        auto file_path = files.create_text_file(content, ".snp.eff");

        REQUIRE_NOTHROW(
            [&]()
            {
                SnpEffectLoader loader(file_path);
                const auto& effects = loader.effects();

                REQUIRE(effects.size() == 1);
                REQUIRE(loader.has_dom_effects() == true);

                auto it = effects.find("rs1301");
                REQUIRE(it != effects.end());
                REQUIRE(it->second.index == 0);
                REQUIRE(it->second.A1 == 'A');
                REQUIRE(it->second.A2 == 'C');
                REQUIRE(it->second.A1freq == 0.25);
                REQUIRE(it->second.add == 0.123);
                REQUIRE(it->second.dom == 0.045);
            }());
    }
}

TEST_CASE("has_dom_effect_column - basic functionality", "[data][snp_effect]")
{
    FileFixture files;

    SECTION("Happy path - file with Dom column")
    {
        std::string content = create_snp_effect_content(
            "ID\tA1\tA2\tA1Frq\tAdd\tDom", {"rs001\tA\tC\t0.25\t0.123\t0.045"});

        auto file_path = files.create_text_file(content, ".snp.eff");

        REQUIRE(has_dom_effect_column(file_path) == true);
    }

    SECTION("Happy path - file without Dom column")
    {
        std::string content = create_snp_effect_content(
            "ID\tA1\tA2\tA1Frq\tAdd", {"rs001\tA\tC\t0.25\t0.123"});

        auto file_path = files.create_text_file(content, ".snp.eff");

        REQUIRE(has_dom_effect_column(file_path) == false);
    }

    SECTION("Happy path - different column order with Dom")
    {
        std::string content = create_snp_effect_content(
            "Dom\tA1Frq\tAdd\tID\tA2\tA1", {"0.045\t0.25\t0.123\trs001\tC\tA"});

        auto file_path = files.create_text_file(content, ".snp.eff");

        REQUIRE(has_dom_effect_column(file_path) == true);
    }

    SECTION("Happy path - file with only header containing Dom")
    {
        std::string content = "ID\tA1\tA2\tA1Frq\tAdd\tDom\n";
        auto file_path = files.create_text_file(content, ".snp.eff");

        REQUIRE(has_dom_effect_column(file_path) == true);
    }

    SECTION("Happy path - file with only header without Dom")
    {
        std::string content = "ID\tA1\tA2\tA1Frq\tAdd\n";
        auto file_path = files.create_text_file(content, ".snp.eff");

        REQUIRE(has_dom_effect_column(file_path) == false);
    }
}
