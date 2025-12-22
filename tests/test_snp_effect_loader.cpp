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
            .chrom = 0,
            .id = 1,
            .pos = 2,
            .a1 = 3,
            .a2 = 4,
            .a1frq = 5,
            .add = 6,
            .dom = 7};

        REQUIRE(indices.has_required_columns() == true);
    }

    SECTION("Happy path - all required columns present, dom optional")
    {
        ColumnIndices indices{
            .chrom = 0,
            .id = 1,
            .pos = 2,
            .a1 = 3,
            .a2 = 4,
            .a1frq = 5,
            .add = 6,
            .dom = -1};

        REQUIRE(indices.has_required_columns() == true);
    }

    SECTION("Exception path - missing ID column")
    {
        ColumnIndices indices{
            .chrom = 0,
            .id = -1,
            .pos = 1,
            .a1 = 2,
            .a2 = 3,
            .a1frq = 4,
            .add = 5,
            .dom = 6};

        REQUIRE(indices.has_required_columns() == false);
    }

    SECTION("Exception path - missing A1 column")
    {
        ColumnIndices indices{
            .chrom = 0,
            .id = 1,
            .pos = 2,
            .a1 = -1,
            .a2 = 3,
            .a1frq = 4,
            .add = 5,
            .dom = 6};

        REQUIRE(indices.has_required_columns() == false);
    }

    SECTION("Exception path - missing A2 column")
    {
        ColumnIndices indices{
            .chrom = 0,
            .id = 1,
            .pos = 2,
            .a1 = 3,
            .a2 = -1,
            .a1frq = 4,
            .add = 5,
            .dom = 6};

        REQUIRE(indices.has_required_columns() == false);
    }

    SECTION("Exception path - missing A1Freq column")
    {
        ColumnIndices indices{
            .chrom = 0,
            .id = 1,
            .pos = 2,
            .a1 = 3,
            .a2 = 4,
            .a1frq = -1,
            .add = 5,
            .dom = 6};

        REQUIRE(indices.has_required_columns() == false);
    }

    SECTION("Exception path - missing Add column")
    {
        ColumnIndices indices{
            .chrom = 0,
            .id = 1,
            .pos = 2,
            .a1 = 3,
            .a2 = 4,
            .a1frq = 5,
            .add = -1,
            .dom = 6};

        REQUIRE(indices.has_required_columns() == false);
    }
}

TEST_CASE("ColumnIndices - max_required_index", "[data][snp_effect]")
{
    SECTION("Happy path - all columns present")
    {
        ColumnIndices indices{
            .chrom = 0,
            .id = 2,
            .pos = 1,
            .a1 = 3,
            .a2 = 4,
            .a1frq = 5,
            .add = 6,
            .dom = 7};

        REQUIRE(indices.max_required_index() == 7);
    }

    SECTION("Happy path - dom column has highest index")
    {
        ColumnIndices indices{
            .chrom = 0,
            .id = 4,
            .pos = 1,
            .a1 = 2,
            .a2 = 3,
            .a1frq = 5,
            .add = 6,
            .dom = 7};

        REQUIRE(indices.max_required_index() == 7);
    }

    SECTION("Happy path - add column has highest index, no dom")
    {
        ColumnIndices indices{
            .chrom = 0,
            .id = 2,
            .pos = 1,
            .a1 = 3,
            .a2 = 4,
            .a1frq = 5,
            .add = 6,
            .dom = -1};

        REQUIRE(indices.max_required_index() == 6);
    }

    SECTION("Happy path - a1frq column has highest index")
    {
        ColumnIndices indices{
            .chrom = 0,
            .id = 2,
            .pos = 1,
            .a1 = 3,
            .a2 = 4,
            .a1frq = 7,
            .add = 5,
            .dom = 6};

        REQUIRE(indices.max_required_index() == 7);
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
            "Chrom\tPosition\tID\tA1\tA2\tA1Freq\tAdd\tDom",
            {"1\t1000\trs001\tA\tC\t0.25\t0.123\t0.045",
             "1\t2000\trs002\tT\tG\t0.75\t-0.456\t0.089",
             "1\t3000\trs003\tC\tA\t0.50\t0.789\t-0.012"});

        auto file_path = files.create_text_file(content, ".snp.eff");

        REQUIRE_NOTHROW(
            [&]()
            {
                SnpEffectLoader loader(file_path);
                const auto& effects = loader.effects();

                REQUIRE(effects.size() == 3);
                REQUIRE(loader.has_dom_effects() == true);

                // Check first SNP
                auto idx1 = effects.find_index("rs001");
                REQUIRE(idx1.has_value());
                REQUIRE(idx1.value() == 0);
                REQUIRE(effects[*idx1].A1 == 'A');
                REQUIRE(effects[*idx1].A2 == 'C');
                REQUIRE(effects[*idx1].chrom == "1");
                REQUIRE(effects[*idx1].pos == 1000);
                REQUIRE(effects.frequencies()[*idx1] == 0.25);
                REQUIRE(effects.additive_effects()[*idx1] == 0.123);
                REQUIRE(effects.dominance_effects()[*idx1] == 0.045);

                // Check second SNP
                auto idx2 = effects.find_index("rs002");
                REQUIRE(idx2.has_value());
                REQUIRE(idx2.value() == 1);
                REQUIRE(effects[*idx2].A1 == 'T');
                REQUIRE(effects[*idx2].A2 == 'G');
                REQUIRE(effects[*idx2].chrom == "1");
                REQUIRE(effects[*idx2].pos == 2000);
                REQUIRE(effects.frequencies()[*idx2] == 0.75);
                REQUIRE(effects.additive_effects()[*idx2] == -0.456);
                REQUIRE(effects.dominance_effects()[*idx2] == 0.089);

                // Check third SNP
                auto idx3 = effects.find_index("rs003");
                REQUIRE(idx3.has_value());
                REQUIRE(idx3.value() == 2);
                REQUIRE(effects[*idx3].A1 == 'C');
                REQUIRE(effects[*idx3].A2 == 'A');
                REQUIRE(effects[*idx3].chrom == "1");
                REQUIRE(effects[*idx3].pos == 3000);
                REQUIRE(effects.frequencies()[*idx3] == 0.50);
                REQUIRE(effects.additive_effects()[*idx3] == 0.789);
                REQUIRE(effects.dominance_effects()[*idx3] == -0.012);
            }());
    }

    SECTION("Happy path - load file without Dom column")
    {
        std::string content = create_snp_effect_content(
            "Chrom\tPosition\tID\tA1\tA2\tA1Freq\tAdd",
            {"1\t1000\trs101\tG\tT\t0.33\t0.111",
             "1\t2000\trs102\tA\tC\t0.67\t-0.222",
             "1\t3000\trs103\tT\tA\t0.90\t0.333"});

        auto file_path = files.create_text_file(content, ".snp.eff");

        REQUIRE_NOTHROW(
            [&]()
            {
                SnpEffectLoader loader(file_path);
                const auto& effects = loader.effects();

                REQUIRE(effects.size() == 3);
                REQUIRE(loader.has_dom_effects() == false);

                // Check first SNP
                auto idx1 = effects.find_index("rs101");
                REQUIRE(idx1.has_value());
                REQUIRE(idx1.value() == 0);
                REQUIRE(effects[*idx1].A1 == 'G');
                REQUIRE(effects[*idx1].A2 == 'T');
                REQUIRE(effects[*idx1].chrom == "1");
                REQUIRE(effects[*idx1].pos == 1000);
                REQUIRE(effects.frequencies()[*idx1] == 0.33);
                REQUIRE(effects.additive_effects()[*idx1] == 0.111);

                auto idx2 = effects.find_index("rs102");
                REQUIRE(idx2.has_value());
                REQUIRE(idx2.value() == 1);
                REQUIRE(effects[*idx2].A1 == 'A');
                REQUIRE(effects[*idx2].A2 == 'C');
                REQUIRE(effects[*idx2].chrom == "1");
                REQUIRE(effects[*idx2].pos == 2000);
                REQUIRE(effects.frequencies()[*idx2] == 0.67);
                REQUIRE(effects.additive_effects()[*idx2] == -0.222);
            }());
    }

    SECTION("Happy path - take_effects moves data")
    {
        std::string content = create_snp_effect_content(
            "Chrom\tPosition\tID\tA1\tA2\tA1Freq\tAdd\tDom",
            {"1\t1000\trs301\tA\tC\t0.25\t0.123\t0.045",
             "1\t2000\trs302\tT\tG\t0.75\t-0.456\t0.089"});

        auto file_path = files.create_text_file(content, ".snp.eff");

        SnpEffectLoader loader(file_path);
        SnpEffects effects = std::move(loader).take_effects();

        REQUIRE(effects.size() == 2);
        REQUIRE(loader.effects().size() == 0);  // Effects should be moved out
    }
}

TEST_CASE("SnpEffectLoader - Error handling", "[data][snp_effect]")
{
    FileFixture files;

    SECTION("Exception path - missing required column in header")
    {
        std::string content = create_snp_effect_content(
            "Chrom\tPosition\tID\tA1\tA2\tA1Freq",  // Missing Add column
            {"1\t1000\trs401\tA\tC\t0.25", "1\t2000\trs402\tT\tG\t0.75"});

        auto file_path = files.create_text_file(content, ".snp.eff");

        REQUIRE_THROWS_MATCHES(
            [&]() { SnpEffectLoader loader(file_path); }(),
            gelex::FileFormatException,
            Catch::Matchers::MessageMatches(EndsWith(
                "missing required columns (ID, Chrom, Pos, A1, A2, A1Freq, "
                "Add)")));
    }

    SECTION("Exception path - insufficient columns in data row")
    {
        std::string content = create_snp_effect_content(
            "Chrom\tPosition\tID\tA1\tA2\tA1Freq\tAdd\tDom",
            {"1\t1000\trs501\tA\tC\t0.25\t0.123\t0.045",
             "1\t2000\trs502\tT\tG\t0.75\t-0.456",  // Missing Dom column
             "1\t3000\trs503\tC\tA\t0.50\t0.789\t-0.012"});

        auto file_path = files.create_text_file(content, ".snp.eff");

        REQUIRE_THROWS_MATCHES(
            [&]() { SnpEffectLoader loader(file_path); }(),
            gelex::FileFormatException,
            Catch::Matchers::MessageMatches(EndsWith(
                "has insufficient columns. Expected at least 8, got 7")));
    }

    SECTION("Exception path - invalid number format in A1Freq")
    {
        std::string content = create_snp_effect_content(
            "Chrom\tPosition\tID\tA1\tA2\tA1Freq\tAdd\tDom",
            {"1\t1000\trs601\tA\tC\t0.25\t0.123\t0.045",
             "1\t2000\trs602\tT\tG\tinvalid\t-0.456\t0.089",  // Invalid A1Freq
             "1\t3000\trs603\tC\tA\t0.50\t0.789\t-0.012"});

        auto file_path = files.create_text_file(content, ".snp.eff");

        REQUIRE_THROWS_MATCHES(
            [&]() { SnpEffectLoader loader(file_path); }(),
            gelex::FileFormatException,
            Catch::Matchers::MessageMatches(EndsWith("as number")));
    }

    SECTION("Exception path - invalid number format in Add")
    {
        std::string content = create_snp_effect_content(
            "Chrom\tPosition\tID\tA1\tA2\tA1Freq\tAdd\tDom",
            {"1\t1000\trs701\tA\tC\t0.25\t0.123\t0.045",
             "1\t2000\trs702\tT\tG\t0.75\tnot_a_number\t0.089",  // Invalid Add
             "1\t3000\trs703\tC\tA\t0.50\t0.789\t-0.012"});

        auto file_path = files.create_text_file(content, ".snp.eff");

        REQUIRE_THROWS_MATCHES(
            [&]() { SnpEffectLoader loader(file_path); }(),
            gelex::FileFormatException,
            Catch::Matchers::MessageMatches(EndsWith("as number")));
    }

    SECTION("Exception path - invalid number format in Dom")
    {
        std::string content = create_snp_effect_content(
            "Chrom\tPosition\tID\tA1\tA2\tA1Freq\tAdd\tDom",
            {"1\t1000\trs801\tA\tC\t0.25\t0.123\t0.045",
             "1\t2000\trs802\tT\tG\t0.75\t-0.456\tinvalid",  // Invalid Dom
             "1\t3000\trs803\tC\tA\t0.50\t0.789\t-0.012"});

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
        std::string content = "Chrom\tPosition\tID\tA1\tA2\tA1Freq\tAdd\tDom\n";
        auto file_path = files.create_text_file(content, ".snp.eff");

        REQUIRE_NOTHROW(
            [&]()
            {
                SnpEffectLoader loader(file_path);
                REQUIRE(loader.effects().size() == 0);
            }());
    }
}

TEST_CASE("SnpEffectLoader - Column order variations", "[data][snp_effect]")
{
    FileFixture files;

    SECTION("Happy path - different column order")
    {
        std::string content = create_snp_effect_content(
            "Chrom\tPosition\tA1Freq\tAdd\tID\tA2\tA1\tDom",  // Different order
            {"1\t1000\t0.25\t0.123\trs1001\tC\tA\t0.045",
             "1\t2000\t0.75\t-0.456\trs1002\tG\tT\t0.089",
             "1\t3000\t0.50\t0.789\trs1003\tA\tC\t-0.012"});

        auto file_path = files.create_text_file(content, ".snp.eff");

        REQUIRE_NOTHROW(
            [&]()
            {
                SnpEffectLoader loader(file_path);
                const auto& effects = loader.effects();

                REQUIRE(effects.size() == 3);
                REQUIRE(loader.has_dom_effects() == true);

                // Check first SNP
                auto idx1 = effects.find_index("rs1001");
                REQUIRE(idx1.has_value());
                REQUIRE(effects[*idx1].A1 == 'A');
                REQUIRE(effects[*idx1].A2 == 'C');
                REQUIRE(effects[*idx1].chrom == "1");
                REQUIRE(effects[*idx1].pos == 1000);
                REQUIRE(effects.frequencies()[*idx1] == 0.25);
                REQUIRE(effects.additive_effects()[*idx1] == 0.123);
                REQUIRE(effects.dominance_effects()[*idx1] == 0.045);
            }());
    }

    SECTION("Happy path - extra columns ignored")
    {
        std::string content = create_snp_effect_content(
            "Chrom\tPosition\tID\tA1\tA2\tA1Freq\tAdd\tExtra3\tDom\tExtra1\tExt"
            "ra2",  // Extra columns
            {"1\t1000\trs1201\tA\tC\t0.25\t0.123\t0.03\t0."
             "045\tignore1\tignore2",
             "1\t2000\trs1202\tT\tG\t0.75\t-0.456\t0.02\t0."
             "089\tignore3\tignore4",
             "1\t3000\trs1203\tC\tA\t0.50\t0.789\t0.03\t-0."
             "012\tignore5\tignore6"});

        auto file_path = files.create_text_file(content, ".snp.eff");

        REQUIRE_NOTHROW(
            [&]()
            {
                SnpEffectLoader loader(file_path);
                const auto& effects = loader.effects();

                REQUIRE(effects.size() == 3);

                auto idx1 = effects.find_index("rs1201");
                REQUIRE(idx1.has_value());
                REQUIRE(effects[*idx1].A1 == 'A');
                REQUIRE(effects[*idx1].A2 == 'C');
                REQUIRE(effects[*idx1].chrom == "1");
                REQUIRE(effects[*idx1].pos == 1000);
                REQUIRE(effects.frequencies()[*idx1] == 0.25);
                REQUIRE(effects.additive_effects()[*idx1] == 0.123);
                REQUIRE(effects.dominance_effects()[*idx1] == 0.045);
            }());
    }
}

TEST_CASE("SnpEffectLoader - Edge cases", "[data][snp_effect]")
{
    FileFixture files;

    SECTION("Happy path - single SNP")
    {
        std::string content = create_snp_effect_content(
            "Chrom\tPosition\tID\tA1\tA2\tA1Freq\tAdd\tDom",
            {"1\t1000\trs1301\tA\tC\t0.25\t0.123\t0.045"});

        auto file_path = files.create_text_file(content, ".snp.eff");

        REQUIRE_NOTHROW(
            [&]()
            {
                SnpEffectLoader loader(file_path);
                const auto& effects = loader.effects();

                REQUIRE(effects.size() == 1);
                REQUIRE(loader.has_dom_effects() == true);

                auto idx = effects.find_index("rs1301");
                REQUIRE(idx.has_value());
                REQUIRE(idx.value() == 0);
                REQUIRE(effects[*idx].A1 == 'A');
                REQUIRE(effects[*idx].A2 == 'C');
                REQUIRE(effects[*idx].chrom == "1");
                REQUIRE(effects[*idx].pos == 1000);
                REQUIRE(effects.frequencies()[*idx] == 0.25);
                REQUIRE(effects.additive_effects()[*idx] == 0.123);
                REQUIRE(effects.dominance_effects()[*idx] == 0.045);
            }());
    }

    SECTION("Happy path - exclude SNPs with nan/inf values")
    {
        std::string content = create_snp_effect_content(
            "Chrom\tPosition\tID\tA1\tA2\tA1Freq\tAdd\tDom",
            {"1\t1000\trs1401\tA\tC\tnan\t0.123\t0.045",     // nan in A1Freq
             "1\t2000\trs1402\tT\tG\t0.75\tInf\t0.089",      // Inf in Add
             "1\t3000\trs1403\tC\tA\t0.50\t0.789\t-Inf",     // -Inf in Dom
             "1\t4000\trs1404\tG\tT\t0.33\t0.111\t0.022"});  // Valid SNP

        auto file_path = files.create_text_file(content, ".snp.eff");

        REQUIRE_NOTHROW(
            [&]()
            {
                SnpEffectLoader loader(file_path);
                const auto& effects = loader.effects();

                // All SNPs are loaded (including NaN/Inf)
                REQUIRE(effects.size() == 1);
                REQUIRE(loader.has_dom_effects() == true);

                auto idx = effects.find_index("rs1404");
                REQUIRE(idx.has_value());
                REQUIRE(idx.value() == 0);  // rs1404 is the 4th SNP, index 3
                REQUIRE(effects[*idx].A1 == 'G');
                REQUIRE(effects[*idx].A2 == 'T');
                REQUIRE(effects[*idx].chrom == "1");
                REQUIRE(effects[*idx].pos == 4000);
                REQUIRE(effects.frequencies()[*idx] == 0.33);
                REQUIRE(effects.additive_effects()[*idx] == 0.111);
                REQUIRE(effects.dominance_effects()[*idx] == 0.022);
            }());
    }
}

TEST_CASE("check_dom_effect_column - basic functionality", "[data][snp_effect]")
{
    FileFixture files;

    SECTION("Happy path - file with Dom column")
    {
        std::string content = create_snp_effect_content(
            "Chrom\tPosition\tID\tA1\tA2\tA1Freq\tAdd\tDom",
            {"1\t1000\trs001\tA\tC\t0.25\t0.123\t0.045"});

        auto file_path = files.create_text_file(content, ".snp.eff");

        REQUIRE(check_dom_effect_column(file_path) == true);
    }

    SECTION("Happy path - file without Dom column")
    {
        std::string content = create_snp_effect_content(
            "ID\tA1\tA2\tA1Freq\tAdd", {"rs001\tA\tC\t0.25\t0.123"});

        auto file_path = files.create_text_file(content, ".snp.eff");

        REQUIRE(check_dom_effect_column(file_path) == false);
    }

    SECTION("Happy path - different column order with Dom")
    {
        std::string content = create_snp_effect_content(
            "Dom\tA1Freq\tAdd\tID\tA2\tA1",
            {"0.045\t0.25\t0.123\trs001\tC\tA"});

        auto file_path = files.create_text_file(content, ".snp.eff");

        REQUIRE(check_dom_effect_column(file_path) == true);
    }

    SECTION("Happy path - file with only header containing Dom")
    {
        std::string content = "ID\tA1\tA2\tA1Freq\tAdd\tDom\n";
        auto file_path = files.create_text_file(content, ".snp.eff");

        REQUIRE(check_dom_effect_column(file_path) == true);
    }

    SECTION("Happy path - file with only header without Dom")
    {
        std::string content = "ID\tA1\tA2\tA1Freq\tAdd\n";
        auto file_path = files.create_text_file(content, ".snp.eff");

        REQUIRE(check_dom_effect_column(file_path) == false);
    }
}
