#include <Eigen/Core>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <memory>
#include <sstream>
#include <string>
#include <string_view>
#include <vector>

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_exception.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>

#include "../src/predict/snp_matcher.h"
#include "bed_fixture.h"
#include "file_fixture.h"
#include "gelex/data/bed_pipe.h"
#include "gelex/data/sample_manager.h"
#include "gelex/exception.h"

namespace fs = std::filesystem;

using namespace gelex;  // NOLINT
using Catch::Matchers::EndsWith;
using Catch::Matchers::WithinAbs;
using gelex::test::are_matrices_equal;
using gelex::test::BedFixture;
using gelex::test::FileFixture;

TEST_CASE("BedPipe - Construction with valid BED files", "[data][bed_pipe]")
{
    BedFixture fixture;

    SECTION("Happy path - dense mapping (all samples match and in same order)")
    {
        auto [bed_prefix, genotypes] = fixture.create_bed_files(10, 20);

        auto fam_path = bed_prefix;
        fam_path.replace_extension(".fam");
        auto sample_manager = std::make_shared<SampleManager>(fam_path, false);
        sample_manager->finalize();

        REQUIRE_NOTHROW(
            [&]()
            {
                BedPipe pipe(bed_prefix, sample_manager);
                REQUIRE(pipe.num_samples() == 10);
                REQUIRE(pipe.num_snps() == 20);
            }());
    }

    SECTION("Happy path - sparse mapping (partial sample overlap)")
    {
        auto [bed_prefix, genotypes] = fixture.create_bed_files(10, 20, 0.1);

        auto fam_path = bed_prefix;
        fam_path.replace_extension(".fam");
        auto sample_manager = std::make_shared<SampleManager>(fam_path, false);

        std::ifstream fam_file(fam_path);
        std::vector<std::string> raw_ids;
        std::string line;
        while (std::getline(fam_file, line))
        {
            std::istringstream iss(line);
            std::string fid;
            std::string iid;
            iss >> fid >> iid;
            raw_ids.push_back(std::format("{}_{}", fid, iid));
        }

        // 只保留前 5 个样本
        std::vector<std::string> intersect_ids;
        for (size_t i = 0; i < 5; ++i)
        {
            intersect_ids.push_back(raw_ids[i]);
        }
        sample_manager->intersect(intersect_ids);
        sample_manager->finalize();

        REQUIRE_NOTHROW(
            [&]()
            {
                BedPipe pipe(bed_prefix, sample_manager);
                REQUIRE(pipe.num_samples() == 5);
                REQUIRE(pipe.num_snps() == 20);
            }());
    }

    SECTION("Exception - file not found")
    {
        BedFixture fixture;
        auto bed_prefix = fixture.create_bed_files(5, 10, 0.0).first;

        auto bed_path = bed_prefix;
        bed_path.replace_extension(".bed");
        fs::remove(bed_path);

        auto fam_path = bed_prefix;
        fam_path.replace_extension(".fam");
        auto sample_manager = std::make_shared<SampleManager>(fam_path, false);
        sample_manager->finalize();

        REQUIRE_THROWS_MATCHES(
            BedPipe(bed_prefix, sample_manager),
            FileOpenException,
            Catch::Matchers::MessageMatches(
                EndsWith("failed to mmap bed file")));
    }

    SECTION("Exception - invalid magic number")
    {
        BedFixture fixture;
        auto bed_prefix = fixture.create_bed_files(5, 10, 0.0).first;

        auto bed_path = bed_prefix;
        bed_path.replace_extension(".bed");
        {
            std::fstream bed_file(
                bed_path, std::ios::in | std::ios::out | std::ios::binary);
            bed_file.seekp(0);
            bed_file.put(0x00);
            bed_file.put(0x00);
            bed_file.put(0x00);
        }

        auto fam_path = bed_prefix;
        fam_path.replace_extension(".fam");
        auto sample_manager = std::make_shared<SampleManager>(fam_path, false);
        sample_manager->finalize();

        REQUIRE_THROWS_MATCHES(
            BedPipe(bed_prefix, sample_manager),
            FileFormatException,
            Catch::Matchers::MessageMatches(
                EndsWith("invalid BED magic number")));
    }

    SECTION("Exception - file too short")
    {
        BedFixture fixture;
        auto bed_prefix = fixture.create_bed_files(5, 10, 0.0).first;

        auto bed_path = bed_prefix;
        bed_path.replace_extension(".bed");
        {
            std::ofstream bed_file(
                bed_path, std::ios::binary | std::ios::trunc);
            bed_file.put(0x6C);
            bed_file.put(0x1B);
            bed_file.put(0x01);
        }

        auto fam_path = bed_prefix;
        fam_path.replace_extension(".fam");
        auto sample_manager = std::make_shared<SampleManager>(fam_path, false);
        sample_manager->finalize();

        REQUIRE_THROWS_AS(
            BedPipe(bed_prefix, sample_manager), FileFormatException);
    }

    SECTION("Exception - sample manager is nullptr")
    {
        BedFixture fixture;
        auto bed_prefix = fixture.create_bed_files(5, 10, 0.0).first;

        REQUIRE_THROWS_MATCHES(
            BedPipe(bed_prefix, nullptr),
            ArgumentValidationException,
            Catch::Matchers::MessageMatches(
                EndsWith("SampleManager cannot be null")));
    }

    SECTION("Exception - BIM file missing")
    {
        BedFixture fixture;
        auto bed_prefix = fixture.create_bed_files(5, 10, 0.0).first;

        auto bim_path = bed_prefix;
        bim_path.replace_extension(".bim");
        fs::remove(bim_path);

        auto fam_path = bed_prefix;
        fam_path.replace_extension(".fam");
        auto sample_manager = std::make_shared<SampleManager>(fam_path, false);
        sample_manager->finalize();

        REQUIRE_THROWS_MATCHES(
            BedPipe(bed_prefix, sample_manager),
            FileNotFoundException,
            Catch::Matchers::MessageMatches(EndsWith("not found")));
    }
}

TEST_CASE("BedPipe - load() method", "[data][bed_pipe]")
{
    BedFixture fixture;

    SECTION("Happy path - load full matrix with dense mapping")
    {
        const Eigen::Index num_samples = 10;
        const Eigen::Index num_snps = 20;
        auto [bed_prefix, genotypes]
            = fixture.create_bed_files(num_samples, num_snps, 0.1);

        auto fam_path = bed_prefix;
        fam_path.replace_extension(".fam");
        auto sample_manager = std::make_shared<SampleManager>(fam_path, true);
        sample_manager->finalize();

        BedPipe pipe(bed_prefix, sample_manager);

        Eigen::MatrixXd loaded = pipe.load();

        REQUIRE(loaded.rows() == num_samples);
        REQUIRE(loaded.cols() == num_snps);
        REQUIRE(are_matrices_equal(loaded, genotypes, 1e-8));
    }

    SECTION("Edge case - minimal samples (1 sample)")
    {
        const Eigen::Index num_samples = 1;
        const Eigen::Index num_snps = 5;
        auto [bed_prefix, genotypes]
            = fixture.create_bed_files(num_samples, num_snps, 0.0);

        auto fam_path = bed_prefix;
        fam_path.replace_extension(".fam");
        auto sample_manager = std::make_shared<SampleManager>(fam_path, false);
        sample_manager->finalize();

        BedPipe pipe(bed_prefix, sample_manager);

        Eigen::MatrixXd loaded = pipe.load();
        REQUIRE(loaded.rows() == num_samples);
        REQUIRE(loaded.cols() == num_snps);
        REQUIRE(are_matrices_equal(loaded, genotypes, 1e-8));
    }

    SECTION("Edge case - minimal SNPs (1 SNP)")
    {
        const Eigen::Index num_samples = 5;
        const Eigen::Index num_snps = 1;
        auto [bed_prefix, genotypes]
            = fixture.create_bed_files(num_samples, num_snps, 0.0);

        auto fam_path = bed_prefix;
        fam_path.replace_extension(".fam");
        auto sample_manager = std::make_shared<SampleManager>(fam_path, false);
        sample_manager->finalize();

        BedPipe pipe(bed_prefix, sample_manager);

        Eigen::MatrixXd loaded = pipe.load();
        REQUIRE(loaded.rows() == num_samples);
        REQUIRE(loaded.cols() == num_snps);
        REQUIRE(are_matrices_equal(loaded, genotypes, 1e-8));
    }

    SECTION("Happy path - sparse mapping with partial samples")
    {
        const Eigen::Index num_raw_samples = 10;
        const Eigen::Index num_snps = 15;
        auto [bed_prefix, genotypes]
            = fixture.create_bed_files(num_raw_samples, num_snps, 0.1);

        auto fam_path = bed_prefix;
        fam_path.replace_extension(".fam");
        auto sample_manager = std::make_shared<SampleManager>(
            fam_path, true);  // set true to avoid reorder

        std::ifstream fam_file(fam_path);
        std::vector<std::string> raw_ids;
        std::string line;
        while (std::getline(fam_file, line))
        {
            std::istringstream iss(line);
            std::string fid;
            std::string iid;
            iss >> fid >> iid;
            raw_ids.push_back(iid);
        }

        std::vector<std::string> intersect_ids(5);
        std::array<Eigen::Index, 5> indices = {0, 2, 4, 6, 8};
        for (auto index : indices)
        {
            intersect_ids.push_back(raw_ids[index]);
        }
        sample_manager->intersect(intersect_ids);
        sample_manager->finalize();

        BedPipe pipe(bed_prefix, sample_manager);

        Eigen::MatrixXd expected(num_raw_samples / 2, num_snps);
        for (Eigen::Index i = 0; i < 5; ++i)
        {
            expected.row(i) = genotypes.row(i * 2);
        }

        Eigen::MatrixXd loaded = pipe.load();
        REQUIRE(loaded.rows() == 5);
        REQUIRE(loaded.cols() == num_snps);
        REQUIRE(are_matrices_equal(loaded, expected, 1e-8));
    }
}

TEST_CASE("BedPipe - load_chunk() method", "[data][bed_pipe]")
{
    BedFixture fixture;

    SECTION("Happy path - load entire matrix as chunk")
    {
        const Eigen::Index num_samples = 8;
        const Eigen::Index num_snps = 12;
        auto [bed_prefix, genotypes]
            = fixture.create_bed_files(num_samples, num_snps, 0.1);

        auto fam_path = bed_prefix;
        fam_path.replace_extension(".fam");
        auto sample_manager = std::make_shared<SampleManager>(fam_path, true);
        sample_manager->finalize();

        BedPipe pipe(bed_prefix, sample_manager);

        Eigen::MatrixXd full_chunk = pipe.load_chunk(0, num_snps);
        Eigen::MatrixXd full_load = pipe.load();

        REQUIRE(full_chunk.rows() == num_samples);
        REQUIRE(full_chunk.cols() == num_snps);
        REQUIRE(full_chunk.rows() == full_load.rows());
        REQUIRE(full_chunk.cols() == full_load.cols());
        REQUIRE(are_matrices_equal(full_chunk, full_load, 1e-8));
        REQUIRE(are_matrices_equal(full_chunk, genotypes, 1e-8));
    }

    SECTION("Happy path - load single column chunk")
    {
        const Eigen::Index num_samples = 6;
        const Eigen::Index num_snps = 10;
        auto [bed_prefix, genotypes]
            = fixture.create_bed_files(num_samples, num_snps, 0.1);

        auto fam_path = bed_prefix;
        fam_path.replace_extension(".fam");
        auto sample_manager = std::make_shared<SampleManager>(fam_path, true);
        sample_manager->finalize();

        BedPipe pipe(bed_prefix, sample_manager);

        Eigen::Index col_idx = 3;
        Eigen::MatrixXd chunk = pipe.load_chunk(col_idx, col_idx + 1);

        REQUIRE(chunk.rows() == num_samples);
        REQUIRE(chunk.cols() == 1);
        REQUIRE(are_matrices_equal(chunk, genotypes.col(col_idx), 1e-8));
    }

    SECTION("Happy path - load middle chunk")
    {
        const Eigen::Index num_samples = 7;
        const Eigen::Index num_snps = 15;
        auto [bed_prefix, genotypes]
            = fixture.create_bed_files(num_samples, num_snps, 0.1);

        auto fam_path = bed_prefix;
        fam_path.replace_extension(".fam");
        auto sample_manager = std::make_shared<SampleManager>(fam_path, true);
        sample_manager->finalize();

        BedPipe pipe(bed_prefix, sample_manager);

        Eigen::Index start = 5;
        Eigen::Index end = 10;
        Eigen::MatrixXd chunk = pipe.load_chunk(start, end);

        REQUIRE(chunk.rows() == num_samples);
        REQUIRE(chunk.cols() == (end - start));

        REQUIRE(are_matrices_equal(
            chunk, genotypes.middleCols(start, end - start), 1e-8));
    }

    SECTION("Exception - invalid chunk range: start < 0")
    {
        const Eigen::Index num_samples = 5;
        const Eigen::Index num_snps = 8;
        auto bed_prefix
            = fixture.create_bed_files(num_samples, num_snps, 0.0).first;

        auto fam_path = bed_prefix;
        fam_path.replace_extension(".fam");
        auto sample_manager = std::make_shared<SampleManager>(fam_path, false);
        sample_manager->finalize();

        BedPipe pipe(bed_prefix, sample_manager);

        REQUIRE_THROWS_AS(pipe.load_chunk(-1, 3), ColumnRangeException);
        REQUIRE_THROWS_AS(
            pipe.load_chunk(0, num_snps + 1), ColumnRangeException);
        REQUIRE_THROWS_AS(
            pipe.load_chunk(3, 3),  // start == end
            ColumnRangeException);
        REQUIRE_THROWS_AS(
            pipe.load_chunk(5, 3),  // start > end
            ColumnRangeException);
    }

    SECTION("Happy path - chunk at beginning")
    {
        const Eigen::Index num_samples = 4;
        const Eigen::Index num_snps = 7;
        auto [bed_prefix, genotypes]
            = fixture.create_bed_files(num_samples, num_snps, 0.1);

        auto fam_path = bed_prefix;
        fam_path.replace_extension(".fam");
        auto sample_manager = std::make_shared<SampleManager>(fam_path, false);
        sample_manager->finalize();

        BedPipe pipe(bed_prefix, sample_manager);

        Eigen::MatrixXd chunk = pipe.load_chunk(0, 3);
        REQUIRE(chunk.rows() == num_samples);
        REQUIRE(chunk.cols() == 3);
        REQUIRE(are_matrices_equal(chunk, genotypes.leftCols(3), 1e-8));
    }

    SECTION("Happy path - chunk at end")
    {
        const Eigen::Index num_samples = 4;
        const Eigen::Index num_snps = 7;
        auto [bed_prefix, genotypes]
            = fixture.create_bed_files(num_samples, num_snps, 0.1);

        auto fam_path = bed_prefix;
        fam_path.replace_extension(".fam");
        auto sample_manager = std::make_shared<SampleManager>(fam_path, false);
        sample_manager->finalize();

        BedPipe pipe(bed_prefix, sample_manager);

        Eigen::MatrixXd chunk = pipe.load_chunk(num_snps - 2, num_snps);
        REQUIRE(chunk.rows() == num_samples);
        REQUIRE(chunk.cols() == 2);
        REQUIRE(are_matrices_equal(chunk, genotypes.rightCols(2), 1e-8));
    }
}

TEST_CASE("BedPipe - sample mapping tests", "[data][bed_pipe]")
{
    BedFixture fixture;

    SECTION("No overlapping samples")
    {
        const Eigen::Index num_raw_samples = 5;
        const Eigen::Index num_snps = 6;
        auto bed_prefix
            = fixture.create_bed_files(num_raw_samples, num_snps, 0.0).first;

        auto fam_path = bed_prefix;
        fam_path.replace_extension(".fam");
        auto sample_manager = std::make_shared<SampleManager>(fam_path, false);

        std::vector<std::string> intersect_ids
            = {"nonexistent_1", "nonexistent_2", "nonexistent_3"};
        sample_manager->intersect(intersect_ids);
        sample_manager->finalize();

        BedPipe pipe(bed_prefix, sample_manager);

        REQUIRE(pipe.num_samples() == 0);

        Eigen::MatrixXd loaded = pipe.load();
        REQUIRE(loaded.rows() == 0);
        REQUIRE(loaded.cols() == num_snps);
    }
}

// TEST_CASE("BedPipe - read plan tests (MatchPlan)", "[data][bed_pipe]")
// {
//     BedFixture fixture;
//
//     SECTION("Happy path - set_read_plan with keep type")
//     {
//         const Eigen::Index num_samples = 5;
//         const Eigen::Index num_snps = 10;
//         auto [bed_prefix, genotypes]
//             = fixture.create_bed_files(num_samples, num_snps, 0.1);
//
//         auto fam_path = bed_prefix;
//         fam_path.replace_extension(".fam");
//         auto sample_manager = std::make_shared<SampleManager>(fam_path,
//         false);
//
//         std::ifstream fam_file(fam_path);
//         std::vector<std::string> raw_ids;
//         std::string line;
//         while (std::getline(fam_file, line))
//         {
//             std::istringstream iss(line);
//             std::string fid, iid;
//             iss >> fid >> iid;
//             raw_ids.push_back(fid + "_" + iid);
//         }
//         std::vector<std::string> intersect_ids;
//         for (size_t i = 0; i < raw_ids.size() - 1; ++i)
//         {
//             intersect_ids.push_back(raw_ids[i]);
//         }
//         sample_manager->intersect(intersect_ids);
//         sample_manager->finalize();
//
//         BedPipe pipe(bed_prefix, sample_manager);
//
//         MatchPlan plan;
//         for (Eigen::Index i = 0; i < num_snps; ++i)
//         {
//             MatchInfo info;
//             if (i == 0 || i == 2 || i == 4)
//             {
//                 info.type = MatchType::keep;
//                 info.target_col = i;
//             }
//             else
//             {
//                 info.type = MatchType::skip;
//                 info.target_col = -1;
//             }
//             plan.push_back(info);
//         }
//
//         pipe.set_read_plan(std::move(plan));
//
//         REQUIRE(pipe.num_snps() == num_snps);
//
//         // 加载数据 - 应该只加载 keep 类型的列
//         Eigen::MatrixXd loaded = pipe.load();
//         REQUIRE(loaded.rows() == pipe.num_samples());
//         REQUIRE(
//             loaded.cols() == num_snps);  // 仍然有 num_snps 列，但跳过的列为
//             NaN
//
//         // 注意：跳过的列可能不是 NaN，如果 is_dense_mapping_ 为 true
//         // 我们只验证 keep 的列有数据
//
//         // 检查 keep 的列是否有数据（可能包含 NaN 由于缺失率）
//         for (Eigen::Index j : {0, 2, 4})
//         {
//             bool has_non_nan = false;
//             for (Eigen::Index i = 0; i < loaded.rows(); ++i)
//             {
//                 if (!std::isnan(loaded(i, j)))
//                 {
//                     has_non_nan = true;
//                     break;
//                 }
//             }
//             REQUIRE(has_non_nan);
//         }
//     }
//
//     SECTION("Happy path - set_read_plan with reverse type")
//     {
//         const Eigen::Index num_samples = 4;
//         const Eigen::Index num_snps = 6;
//         auto bed_prefix = fixture.create_bed_files(num_samples, num_snps,
//         0.0)
//                               .first;  // 无缺失
//
//         auto fam_path = bed_prefix;
//         fam_path.replace_extension(".fam");
//         auto sample_manager = std::make_shared<SampleManager>(fam_path,
//         false); sample_manager->finalize();
//
//         BedPipe pipe(bed_prefix, sample_manager);
//
//         // 首先加载原始数据作为参考
//         Eigen::MatrixXd original = pipe.load();
//
//         // 创建反转第 1 列的计划
//         MatchPlan plan;
//         for (Eigen::Index i = 0; i < num_snps; ++i)
//         {
//             MatchInfo info;
//             if (i == 1)
//             {
//                 info.type = MatchType::reverse;
//                 info.target_col = i;
//             }
//             else
//             {
//                 info.type = MatchType::keep;
//                 info.target_col = i;
//             }
//             plan.push_back(info);
//         }
//
//         pipe.set_read_plan(std::move(plan));
//
//         // 加载应用计划后的数据
//         Eigen::MatrixXd with_plan = pipe.load();
//
//         // 验证维度
//         REQUIRE(with_plan.rows() == num_samples);
//         REQUIRE(with_plan.cols() == num_snps);
//
//         // 检查第 1 列是否反转
//         // 反转映射：0.0 <-> 2.0, 1.0 保持不变
//         for (Eigen::Index i = 0; i < num_samples; ++i)
//         {
//             double orig = original(i, 1);
//             double planned = with_plan(i, 1);
//
//             if (std::abs(orig - 0.0) < 1e-10)
//             {
//                 REQUIRE(std::abs(planned - 2.0) < 1e-10);
//             }
//             else if (std::abs(orig - 2.0) < 1e-10)
//             {
//                 REQUIRE(std::abs(planned - 0.0) < 1e-10);
//             }
//             else if (std::abs(orig - 1.0) < 1e-10)
//             {
//                 REQUIRE(std::abs(planned - 1.0) < 1e-10);
//             }
//         }
//
//         // 检查其他列保持不变
//         for (Eigen::Index j : {0, 2, 3, 4, 5})
//         {
//             for (Eigen::Index i = 0; i < num_samples; ++i)
//             {
//                 REQUIRE(original(i, j) == with_plan(i, j));
//             }
//         }
//     }
//
//     SECTION("Happy path - reset_to_default clears plan")
//     {
//         const Eigen::Index num_samples = 3;
//         const Eigen::Index num_snps = 5;
//         auto bed_prefix
//             = fixture.create_bed_files(num_samples, num_snps, 0.0).first;
//
//         auto fam_path = bed_prefix;
//         fam_path.replace_extension(".fam");
//         auto sample_manager = std::make_shared<SampleManager>(fam_path,
//         false); sample_manager->finalize();
//
//         BedPipe pipe(bed_prefix, sample_manager);
//
//         // 设置计划
//         MatchPlan plan(num_snps);
//         for (Eigen::Index i = 0; i < num_snps; ++i)
//         {
//             plan[i].type = MatchType::skip;
//             plan[i].target_col = -1;
//         }
//         pipe.set_read_plan(std::move(plan));
//
//         // 验证计划已设置
//         Eigen::MatrixXd with_plan = pipe.load();
//         // 维度应该正确
//         REQUIRE(with_plan.rows() == num_samples);
//         REQUIRE(with_plan.cols() == num_snps);
//         // 注意：当 is_dense_mapping_ 为 true 时，跳过的列可能不会被初始化为
//         NaN
//
//         // 重置为默认
//         pipe.reset_to_default();
//
//         // 现在应该加载实际数据
//         Eigen::MatrixXd after_reset = pipe.load();
//         bool has_non_nan = false;
//         for (Eigen::Index j = 0; j < num_snps; ++j)
//         {
//             for (Eigen::Index i = 0; i < num_samples; ++i)
//             {
//                 if (!std::isnan(after_reset(i, j)))
//                 {
//                     has_non_nan = true;
//                     break;
//                 }
//             }
//             if (has_non_nan)
//                 break;
//         }
//         REQUIRE(has_non_nan);
//     }
//
//     SECTION("Happy path - plan with reordering")
//     {
//         const Eigen::Index num_samples = 4;
//         const Eigen::Index num_snps = 8;
//         auto bed_prefix
//             = fixture.create_bed_files(num_samples, num_snps, 0.0).first;
//
//         auto fam_path = bed_prefix;
//         fam_path.replace_extension(".fam");
//         auto sample_manager = std::make_shared<SampleManager>(fam_path,
//         false); sample_manager->finalize();
//
//         BedPipe pipe(bed_prefix, sample_manager);
//
//         // 首先加载原始数据作为参考
//         Eigen::MatrixXd original = pipe.load();
//
//         // 创建重新排序的计划：将列 3, 1, 5 映射到位置 0, 1, 2
//         MatchPlan plan(num_snps);
//         for (Eigen::Index i = 0; i < num_snps; ++i)
//         {
//             plan[i].type = MatchType::skip;
//             plan[i].target_col = -1;
//         }
//
//         // 设置特定的映射
//         plan[3].type = MatchType::keep;
//         plan[3].target_col = 0;
//
//         plan[1].type = MatchType::keep;
//         plan[1].target_col = 1;
//
//         plan[5].type = MatchType::keep;
//         plan[5].target_col = 2;
//
//         pipe.set_read_plan(std::move(plan));
//
//         // 加载数据
//         Eigen::MatrixXd with_plan = pipe.load();
//
//         // 验证维度
//         REQUIRE(with_plan.rows() == num_samples);
//         REQUIRE(with_plan.cols() == num_snps);
//
//         // 检查映射
//         // 列 3 应该出现在位置 0
//         for (Eigen::Index i = 0; i < num_samples; ++i)
//         {
//             REQUIRE(with_plan(i, 0) == original(i, 3));
//         }
//
//         // 列 1 应该出现在位置 1
//         for (Eigen::Index i = 0; i < num_samples; ++i)
//         {
//             REQUIRE(with_plan(i, 1) == original(i, 1));
//         }
//
//         // 列 5 应该出现在位置 2
//         for (Eigen::Index i = 0; i < num_samples; ++i)
//         {
//             REQUIRE(with_plan(i, 2) == original(i, 5));
//         }
//
//         // 注意：跳过的列可能不是 NaN，如果 is_dense_mapping_ 为 true
//         // 我们只验证了映射的列正确
//     }
// }
//
// TEST_CASE("BedPipe - reverse decoding tests", "[data][bed_pipe]")
// {
//     // 注意：反向解码已经在 read plan 测试中部分测试过
//     // 这里进行更全面的测试
//
//     SECTION("Verify reverse LUT mapping")
//     {
//         // 验证反向查找表的映射关系
//         // 标准映射：00->2.0, 01->NaN, 10->1.0, 11->0.0
//         // 反向映射：00->0.0, 01->NaN, 10->1.0, 11->2.0
//
//         // 由于 LUT 是私有的，我们通过加载实际数据来测试
//         BedFixture fixture;
//         const Eigen::Index num_samples = 4;
//         const Eigen::Index num_snps = 1;
//
//         // 创建一个简单的测试矩阵，包含所有可能的基因型
//         Eigen::MatrixXd test_matrix(num_samples, 1);
//         test_matrix << 0.0, 1.0, 2.0,
//         std::numeric_limits<double>::quiet_NaN();
//
//         // 使用 BedFixture 从矩阵创建 BED 文件
//         auto bed_prefix = fixture.create_bed_files_from_matrix(test_matrix);
//
//         auto fam_path = bed_prefix;
//         fam_path.replace_extension(".fam");
//         auto sample_manager = std::make_shared<SampleManager>(fam_path,
//         false); sample_manager->finalize();
//
//         BedPipe pipe(bed_prefix, sample_manager);
//
//         // 加载原始数据
//         Eigen::MatrixXd original = pipe.load();
//
//         // 验证原始数据（应该与输入匹配，除了 NaN）
//         REQUIRE(std::abs(original(0, 0) - 0.0) < 1e-10);
//         REQUIRE(std::abs(original(1, 0) - 1.0) < 1e-10);
//         REQUIRE(std::abs(original(2, 0) - 2.0) < 1e-10);
//         REQUIRE(std::isnan(original(3, 0)));
//
//         // 创建反转计划
//         MatchPlan plan(1);
//         plan[0].type = MatchType::reverse;
//         plan[0].target_col = 0;
//
//         pipe.set_read_plan(std::move(plan));
//
//         // 加载反转后的数据
//         Eigen::MatrixXd reversed = pipe.load();
//
//         // 验证反转映射
//         // 0.0 -> 2.0
//         REQUIRE(std::abs(reversed(0, 0) - 2.0) < 1e-10);
//         // 1.0 -> 1.0 (不变)
//         REQUIRE(std::abs(reversed(1, 0) - 1.0) < 1e-10);
//         // 2.0 -> 0.0
//         REQUIRE(std::abs(reversed(2, 0) - 0.0) < 1e-10);
//         // NaN -> NaN
//         REQUIRE(std::isnan(reversed(3, 0)));
//     }
// }
//
// TEST_CASE("BedPipe - static format_bed_path() method", "[data][bed_pipe]")
// {
//     SECTION("Happy path - path with .bed extension")
//     {
//         BedFixture fixture;
//         auto bed_prefix = fixture.create_bed_files(3, 4, 0.0).first;
//         auto bed_path = bed_prefix;
//         bed_path.replace_extension(".bed");
//
//         // 应该返回相同的路径
//         auto formatted = BedPipe::format_bed_path(bed_path.string());
//         REQUIRE(formatted == bed_path);
//     }
//
//     SECTION("Happy path - path without .bed extension")
//     {
//         BedFixture fixture;
//         auto bed_prefix = fixture.create_bed_files(3, 4, 0.0).first;
//
//         // 传入不带扩展名的路径
//         auto formatted = BedPipe::format_bed_path(bed_prefix.string());
//         REQUIRE(formatted == bed_prefix.replace_extension(".bed"));
//     }
//
//     SECTION("Exception - file not found")
//     {
//         REQUIRE_THROWS_MATCHES(
//             BedPipe::format_bed_path("non_existent_file"),
//             FileNotFoundException,
//             Catch::Matchers::MessageMatches(EndsWith("file not found")));
//
//         REQUIRE_THROWS_MATCHES(
//             BedPipe::format_bed_path("non_existent_file.bed"),
//             FileNotFoundException,
//             Catch::Matchers::MessageMatches(EndsWith("file not found")));
//     }
// }
//
// TEST_CASE("BedPipe - num_samples() and num_snps() methods",
// "[data][bed_pipe]")
// {
//     BedFixture fixture;
//
//     SECTION("Happy path - basic counts")
//     {
//         const Eigen::Index num_samples = 7;
//         const Eigen::Index num_snps = 9;
//         auto bed_prefix
//             = fixture.create_bed_files(num_samples, num_snps, 0.1).first;
//
//         auto fam_path = bed_prefix;
//         fam_path.replace_extension(".fam");
//         auto sample_manager = std::make_shared<SampleManager>(fam_path,
//         false); sample_manager->finalize();
//
//         BedPipe pipe(bed_prefix, sample_manager);
//
//         REQUIRE(pipe.num_samples() == num_samples);
//         REQUIRE(pipe.num_snps() == num_snps);
//     }
//
//     SECTION("num_snps() with custom plan")
//     {
//         const Eigen::Index num_samples = 5;
//         const Eigen::Index num_snps = 6;
//         auto bed_prefix
//             = fixture.create_bed_files(num_samples, num_snps, 0.0).first;
//
//         auto fam_path = bed_prefix;
//         fam_path.replace_extension(".fam");
//         auto sample_manager = std::make_shared<SampleManager>(fam_path,
//         false); sample_manager->finalize();
//
//         BedPipe pipe(bed_prefix, sample_manager);
//
//         // 默认情况下，num_snps() 返回原始 SNP 数
//         REQUIRE(pipe.num_snps() == num_snps);
//
//         // 设置计划
//         MatchPlan plan(num_snps);
//         for (Eigen::Index i = 0; i < num_snps; ++i)
//         {
//             plan[i].type = MatchType::keep;
//             plan[i].target_col = i;
//         }
//         pipe.set_read_plan(std::move(plan));
//
//         // 设置计划后，num_snps() 应该返回计划的大小
//         REQUIRE(pipe.num_snps() == num_snps);
//
//         // 重置后应该恢复
//         pipe.reset_to_default();
//         REQUIRE(pipe.num_snps() == num_snps);
//     }
//
//     SECTION("num_samples() with sparse mapping")
//     {
//         const Eigen::Index num_raw_samples = 10;
//         const Eigen::Index num_snps = 5;
//         auto bed_prefix
//             = fixture.create_bed_files(num_raw_samples, num_snps, 0.1).first;
//
//         // 使用原始 FAM 文件创建 SampleManager
//         auto fam_path = bed_prefix;
//         fam_path.replace_extension(".fam");
//         auto sample_manager = std::make_shared<SampleManager>(fam_path,
//         false);
//
//         // 读取原始样本 ID 用于交集
//         std::ifstream fam_file(fam_path);
//         std::vector<std::string> raw_ids;
//         std::string line;
//         while (std::getline(fam_file, line))
//         {
//             std::istringstream iss(line);
//             std::string fid, iid;
//             iss >> fid >> iid;
//             raw_ids.push_back(fid + "_" + iid);
//         }
//
//         // 只保留前 3 个样本
//         std::vector<std::string> intersect_ids;
//         for (size_t i = 0; i < 3; ++i)
//         {
//             intersect_ids.push_back(raw_ids[i]);
//         }
//         sample_manager->intersect(intersect_ids);
//         sample_manager->finalize();
//
//         BedPipe pipe(bed_prefix, sample_manager);
//
//         // num_samples() 应该返回 sample_manager 中的样本数，而不是原始样本数
//         REQUIRE(pipe.num_samples() == 3);
//     }
// }
//
// TEST_CASE("BedPipe - integration tests", "[data][bed_pipe]")
// {
//     BedFixture fixture;
//
//     SECTION("Complete workflow: load chunk with custom plan")
//     {
//         const Eigen::Index num_samples = 6;
//         const Eigen::Index num_snps = 12;
//         auto bed_prefix
//             = fixture.create_bed_files(num_samples, num_snps, 0.1).first;
//
//         auto fam_path = bed_prefix;
//         fam_path.replace_extension(".fam");
//         auto sample_manager = std::make_shared<SampleManager>(fam_path,
//         false); sample_manager->finalize();
//
//         BedPipe pipe(bed_prefix, sample_manager);
//
//         // 创建计划：只处理偶数索引的 SNP，并反转其中的一半
//         MatchPlan plan(num_snps);
//         for (Eigen::Index i = 0; i < num_snps; ++i)
//         {
//             if (i % 2 == 0)
//             {  // 偶数索引
//                 if (i % 4 == 0)
//                 {  // 每两个偶数中的一个反转
//                     plan[i].type = MatchType::reverse;
//                 }
//                 else
//                 {
//                     plan[i].type = MatchType::keep;
//                 }
//                 plan[i].target_col = i;
//             }
//             else
//             {
//                 plan[i].type = MatchType::skip;
//                 plan[i].target_col = -1;
//             }
//         }
//
//         pipe.set_read_plan(std::move(plan));
//
//         // 加载完整矩阵
//         Eigen::MatrixXd full_with_plan = pipe.load();
//
//         // 注意：load_chunk() 与自定义计划的集成测试暂时跳过
//         // 因为当前实现中 target_col 映射在分块加载时可能需要调整
//     }
//
//     SECTION("Multiple operations on same pipe")
//     {
//         const Eigen::Index num_samples = 4;
//         const Eigen::Index num_snps = 7;
//         auto bed_prefix
//             = fixture.create_bed_files(num_samples, num_snps, 0.0).first;
//
//         auto fam_path = bed_prefix;
//         fam_path.replace_extension(".fam");
//         auto sample_manager = std::make_shared<SampleManager>(fam_path,
//         false); sample_manager->finalize();
//
//         BedPipe pipe(bed_prefix, sample_manager);
//
//         // 操作 1: 加载完整矩阵
//         Eigen::MatrixXd m1 = pipe.load();
//
//         // 操作 2: 加载分块
//         Eigen::MatrixXd m2 = pipe.load_chunk(1, 4);
//
//         // 操作 3: 设置计划并加载
//         MatchPlan plan(num_snps);
//         for (Eigen::Index i = 0; i < num_snps; ++i)
//         {
//             plan[i].type = MatchType::keep;
//             plan[i].target_col = i;
//         }
//         plan[2].type = MatchType::reverse;
//         pipe.set_read_plan(std::move(plan));
//         Eigen::MatrixXd m3 = pipe.load();
//
//         // 操作 4: 重置并再次加载
//         pipe.reset_to_default();
//         Eigen::MatrixXd m4 = pipe.load();
//
//         // 验证 m1 和 m4 应该相同
//         REQUIRE(m1.rows() == m4.rows());
//         REQUIRE(m1.cols() == m4.cols());
//         for (Eigen::Index j = 0; j < num_snps; ++j)
//         {
//             for (Eigen::Index i = 0; i < num_samples; ++i)
//             {
//                 REQUIRE(m1(i, j) == m4(i, j));
//             }
//         }
//     }
// }
