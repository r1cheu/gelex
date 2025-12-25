#include <cmath>
#include <fstream>
#include <map>
#include <sstream>
#include <string>
#include <string_view>
#include <vector>

#include <Eigen/Core>
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_exception.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>

#include "../src/predict/predict_engine.h"
#include "bed_fixture.h"
#include "file_fixture.h"
#include "gelex/exception.h"

namespace fs = std::filesystem;

using namespace gelex;  // NOLINT
using Catch::Matchers::EndsWith;
using Catch::Matchers::WithinAbs;
using gelex::test::BedFixture;
using gelex::test::FileFixture;

namespace
{

std::pair<std::vector<std::string>, std::vector<std::string>> read_fam(
    const std::filesystem::path& path)
{
    std::vector<std::string> fids;
    std::vector<std::string> iids;
    {
        std::ifstream fam_file(path);
        std::string line;
        while (std::getline(fam_file, line))
        {
            std::istringstream iss(line);
            std::string buffer;
            iss >> buffer;  // FID
            fids.push_back(buffer);
            iss >> buffer;  // IID
            iids.push_back(buffer);
        }
    }
    return {std::move(fids), std::move(iids)};
}

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

}  // namespace

TEST_CASE(
    "PredictEngine - Basic prediction with BED + SNP effects",
    "[predict][predict_engine]")
{
    BedFixture bed_fixture;
    FileFixture file_fixture;

    SECTION("Happy path - additive effects only")
    {
        const Eigen::Index num_samples = 5;
        const Eigen::Index num_snps = 10;
        auto [bed_prefix, expected_genotypes]
            = bed_fixture.create_bed_files(num_samples, num_snps, 0.0);

        auto fam_path = bed_prefix;
        fam_path.replace_extension(".fam");
        auto [fids, iids] = read_fam(fam_path);

        std::string snp_header = "Chrom\tPosition\tID\tA1\tA2\tA1Freq\tAdd";
        std::vector<std::string> snp_rows;
        for (Eigen::Index i = 0; i < num_snps; ++i)
        {
            double freq = 0.3 + 0.1 * (i % 5);
            double effect = 0.01 * (i % 3);
            snp_rows.push_back(
                std::format(
                    "1\t{}\trs{}\tA\tC\t{:.6f}\t{:.6f}",
                    1000 + i * 100,
                    i + 1,
                    freq,
                    effect));
        }
        auto snp_path = file_fixture.create_text_file(
            create_snp_effect_content(snp_header, snp_rows), ".snp.eff");

        // 创建协变量效应文件（仅截距）
        auto covar_path = file_fixture.create_text_file(
            "term\tmean\tstddev\tpercentile_5\tpercentile_95\tess\trhat\n"
            "Intercept\t1.5\t0.1\t1.3\t1.7\t1000\t1.0\n",
            ".param");

        // 创建输出文件路径
        auto output_path
            = file_fixture.generate_random_file_path(".predictions");

        // 配置并运行PredictEngine
        PredictEngine::Config config;
        config.bed_path = bed_prefix;
        config.snp_effect_path = snp_path;
        config.covar_effect_path = covar_path;
        config.qcovar_path = "";
        config.dcovar_path = "";
        config.output_path = output_path;
        config.iid_only = false;

        PredictEngine engine(config);

        REQUIRE_NOTHROW(engine.run());

        // 验证预测结果
        const auto& predictions = engine.predictions();
        const auto& sample_ids = engine.sample_ids();

        REQUIRE(predictions.size() == num_samples);
        REQUIRE(sample_ids.size() == static_cast<size_t>(num_samples));

        // 验证样本ID格式
        for (size_t i = 0; i < sample_ids.size(); ++i)
        {
            REQUIRE(sample_ids[i] == std::format("{}_{}", fids[i], iids[i]));
        }

        // 验证输出文件存在
        REQUIRE(fs::exists(output_path));

        // 可选：读取输出文件验证格式
        std::ifstream out_file(output_path);
        std::string header_line;
        std::getline(out_file, header_line);
        REQUIRE(header_line.find("sample_id") != std::string::npos);
        REQUIRE(header_line.find("total_prediction") != std::string::npos);
        REQUIRE(header_line.find("Intercept") != std::string::npos);
        REQUIRE(header_line.find("snp_contribution") != std::string::npos);
    }
}

TEST_CASE(
    "PredictEngine - Full prediction with covariates",
    "[predict][predict_engine]")
{
    BedFixture bed_fixture;
    FileFixture file_fixture;

    SECTION("Happy path with quantitative covariates")
    {
        const Eigen::Index num_samples = 3;
        const Eigen::Index num_snps = 5;
        auto [bed_prefix, expected_genotypes]
            = bed_fixture.create_bed_files(num_samples, num_snps, 0.0);

        // 读取FAM文件
        auto fam_path = bed_prefix;
        fam_path.replace_extension(".fam");
        auto [fids, iids] = read_fam(fam_path);

        // 创建SNP效应文件
        std::string snp_header = "Chrom\tPosition\tID\tA1\tA2\tA1Freq\tAdd";
        std::vector<std::string> snp_rows;
        for (Eigen::Index i = 0; i < num_snps; ++i)
        {
            double freq = 0.25;
            double effect = 0.02 * (i % 2 == 0 ? 1 : -1);
            snp_rows.push_back(
                std::format(
                    "1\t{}\trs{}\tA\tC\t{:.6f}\t{:.6f}",
                    1000 + i * 100,
                    i + 1,
                    freq,
                    effect));
        }
        auto snp_path = file_fixture.create_text_file(
            create_snp_effect_content(snp_header, snp_rows), ".snp.eff");

        // 创建定量协变量文件
        std::string qcovar_content = "FID\tIID\tAge\n";
        for (size_t i = 0; i < fids.size(); ++i)
        {
            qcovar_content
                += std::format("{}\t{}\t{}\n", fids[i], iids[i], 20 + i);
        }
        auto qcovar_path
            = file_fixture.create_text_file(qcovar_content, ".qcovar");

        // 创建协变量效应文件（截距 + Age）
        auto covar_path = file_fixture.create_text_file(
            "term\tmean\tstddev\tpercentile_5\tpercentile_95\tess\trhat\n"
            "Intercept\t1.0\t0.1\t0.8\t1.2\t1000\t1.0\n"
            "Age\t0.3\t0.05\t0.2\t0.4\t800\t1.01\n",
            ".param");

        // 输出文件路径
        auto output_path
            = file_fixture.generate_random_file_path(".predictions");

        PredictEngine::Config config;
        config.bed_path = bed_prefix;
        config.snp_effect_path = snp_path;
        config.covar_effect_path = covar_path;
        config.qcovar_path = qcovar_path;
        config.dcovar_path = "";
        config.output_path = output_path;
        config.iid_only = false;

        PredictEngine engine(config);

        REQUIRE_NOTHROW(engine.run());

        const auto& predictions = engine.predictions();
        const auto& sample_ids = engine.sample_ids();

        REQUIRE(predictions.size() == num_samples);
        REQUIRE(sample_ids.size() == static_cast<size_t>(num_samples));

        // 验证输出文件
        REQUIRE(fs::exists(output_path));

        // 读取输出文件验证协变量列
        std::ifstream out_file(output_path);
        std::string header_line;
        std::getline(out_file, header_line);
        REQUIRE(header_line.find("Age") != std::string::npos);
    }
}

TEST_CASE("PredictEngine - Error handling", "[predict][predict_engine]")
{
    BedFixture bed_fixture;
    FileFixture file_fixture;

    SECTION("Missing output path")
    {
        const Eigen::Index num_samples = 2;
        const Eigen::Index num_snps = 3;
        auto [bed_prefix, _]
            = bed_fixture.create_bed_files(num_samples, num_snps, 0.0);

        auto snp_path = file_fixture.create_text_file(
            "Chrom\tPosition\tID\tA1\tA2\tA1Freq\tAdd\n"
            "1\t1000\trs1\tA\tC\t0.3\t0.1\n",
            ".snp.eff");

        auto covar_path = file_fixture.create_text_file(
            "term\tmean\tstddev\tpercentile_5\tpercentile_95\tess\trhat\n"
            "Intercept\t1.0\t0.1\t0.8\t1.2\t1000\t1.0\n",
            ".param");

        PredictEngine::Config config;
        config.bed_path = bed_prefix;
        config.snp_effect_path = snp_path;
        config.covar_effect_path = covar_path;
        config.qcovar_path = "";
        config.dcovar_path = "";
        config.output_path = "";  // 缺失输出路径
        config.iid_only = false;

        REQUIRE_THROWS_AS(PredictEngine(config), InvalidInputException);
    }

    SECTION("Missing required covariate coefficient")
    {
        const Eigen::Index num_samples = 2;
        const Eigen::Index num_snps = 3;
        auto [bed_prefix, _]
            = bed_fixture.create_bed_files(num_samples, num_snps, 0.0);

        auto snp_path = file_fixture.create_text_file(
            "Chrom\tPosition\tID\tA1\tA2\tA1Freq\tAdd\n"
            "1\t1000\trs1\tA\tC\t0.3\t0.1\n"
            "1\t2000\trs2\tT\tG\t0.4\t-0.2\n"
            "1\t3000\trs3\tC\tA\t0.5\t0.05\n",
            ".snp.eff");

        // 创建定量协变量文件
        std::string qcovar_content
            = "FID\tIID\tAge\n"
              "fid1\tiid1\t30\n"
              "fid1\tiid2\t35\n";
        auto qcovar_path
            = file_fixture.create_text_file(qcovar_content, ".qcovar");

        // 协变量效应文件缺少Age系数
        auto covar_path = file_fixture.create_text_file(
            "term\tmean\tstddev\tpercentile_5\tpercentile_95\tess\trhat\n"
            "Intercept\t1.0\t0.1\t0.8\t1.2\t1000\t1.0\n",
            ".param");

        auto output_path
            = file_fixture.generate_random_file_path(".predictions");

        PredictEngine::Config config;
        config.bed_path = bed_prefix;
        config.snp_effect_path = snp_path;
        config.covar_effect_path = covar_path;
        config.qcovar_path = qcovar_path;
        config.dcovar_path = "";
        config.output_path = output_path;
        config.iid_only = false;

        PredictEngine engine(config);
        REQUIRE_THROWS_AS(engine.run(), InvalidInputException);
    }
}

TEST_CASE("PredictEngine - iid_only mode", "[predict][predict_engine]")
{
    BedFixture bed_fixture;
    FileFixture file_fixture;

    SECTION("iid_only = true")
    {
        const Eigen::Index num_samples = 3;
        const Eigen::Index num_snps = 5;
        auto [bed_prefix, expected_genotypes]
            = bed_fixture.create_bed_files(num_samples, num_snps, 0.0);

        // 读取FAM文件获取IID
        auto fam_path = bed_prefix;
        fam_path.replace_extension(".fam");
        auto [_, iids] = read_fam(fam_path);

        // 创建SNP效应文件
        std::string snp_header = "Chrom\tPosition\tID\tA1\tA2\tA1Freq\tAdd";
        std::vector<std::string> snp_rows;
        for (Eigen::Index i = 0; i < num_snps; ++i)
        {
            snp_rows.push_back(
                std::format(
                    "1\t{}\trs{}\tA\tC\t0.3\t0.01", 1000 + i * 100, i + 1));
        }
        auto snp_path = file_fixture.create_text_file(
            create_snp_effect_content(snp_header, snp_rows), ".snp.eff");

        // 创建定量协变量文件（仅使用IID）
        std::string qcovar_content = "FID\tIID\tAge\n";
        for (size_t i = 0; i < iids.size(); ++i)
        {
            qcovar_content += std::format("1\t{}\t{}\n", iids[i], 20 + i);
        }
        auto qcovar_path
            = file_fixture.create_text_file(qcovar_content, ".qcovar");

        // 协变量效应文件
        auto covar_path = file_fixture.create_text_file(
            "term\tmean\tstddev\tpercentile_5\tpercentile_95\tess\trhat\n"
            "Intercept\t1.0\t0.1\t0.8\t1.2\t1000\t1.0\n"
            "Age\t0.2\t0.05\t0.1\t0.3\t800\t1.01\n",
            ".param");

        auto output_path
            = file_fixture.generate_random_file_path(".predictions");

        PredictEngine::Config config;
        config.bed_path = bed_prefix;
        config.snp_effect_path = snp_path;
        config.covar_effect_path = covar_path;
        config.qcovar_path = qcovar_path;
        config.dcovar_path = "";
        config.output_path = output_path;
        config.iid_only = true;  // 使用IID-only模式

        PredictEngine engine(config);

        REQUIRE_NOTHROW(engine.run());

        const auto& sample_ids = engine.sample_ids();
        REQUIRE(sample_ids.size() == static_cast<size_t>(num_samples));

        // 验证样本ID仅为IID
        for (size_t i = 0; i < sample_ids.size(); ++i)
        {
            REQUIRE(sample_ids[i] == iids[i]);
        }
    }
}

TEST_CASE(
    "PredictEngine - Numerical correctness verification",
    "[predict][predict_engine]")
{
    BedFixture bed_fixture;
    FileFixture file_fixture;

    SECTION("Additive effects only - 2 samples x 2 SNPs")
    {
        // 确定性基因型矩阵
        Eigen::MatrixXd genotypes(2, 2);
        genotypes << 0.0, 2.0, 1.0, 1.0;

        // 创建BED文件
        auto [bed_prefix, _] = bed_fixture.create_deterministic_bed_files(
            genotypes,
            {"sample1", "sample2"},   // 样本ID
            {"rs1", "rs2"},           // SNP ID
            {"1", "1"},               // 染色体
            {{'A', 'C'}, {'T', 'G'}}  // 等位基因对
        );

        // SNP效应文件：频率 p1=0.3, p2=0.4; 加性效应 β1=0.1, β2=-0.05
        std::string snp_header = "Chrom\tPosition\tID\tA1\tA2\tA1Freq\tAdd";
        std::vector<std::string> snp_rows
            = {"1\t1000\trs1\tA\tC\t0.300000\t0.100000",
               "1\t2000\trs2\tT\tG\t0.400000\t-0.050000"};
        auto snp_path = file_fixture.create_text_file(
            create_snp_effect_content(snp_header, snp_rows), ".snp.eff");

        // 协变量效应文件：仅截距 1.0
        auto covar_path = file_fixture.create_text_file(
            "term\tmean\tstddev\tpercentile_5\tpercentile_95\tess\trhat\n"
            "Intercept\t1.0\t0.1\t0.8\t1.2\t1000\t1.0\n",
            ".param");

        auto output_path
            = file_fixture.generate_random_file_path(".predictions");

        PredictEngine::Config config;
        config.bed_path = bed_prefix;
        config.snp_effect_path = snp_path;
        config.covar_effect_path = covar_path;
        config.qcovar_path = "";
        config.dcovar_path = "";
        config.output_path = output_path;
        config.iid_only = false;

        PredictEngine engine(config);
        REQUIRE_NOTHROW(engine.run());

        const auto& predictions = engine.predictions();
        const auto& snp_contributions = engine.snp_predictions();
        const auto& covar_predictions = engine.covar_predictions();

        REQUIRE(predictions.size() == 2);
        REQUIRE(snp_contributions.size() == 2);
        REQUIRE(covar_predictions.rows() == 2);
        REQUIRE(covar_predictions.cols() == 1);  // 仅截距

        // 精确计算期望值（使用与 PredictEngine 相同的公式）
        const double intercept = 1.0;
        const double eps = 1e-10;

        // SNP 参数
        const std::vector<double> frequencies = {0.3, 0.4};
        const std::vector<double> effects = {0.1, -0.05};

        // 基因型矩阵
        const Eigen::MatrixXd& geno = genotypes;

        // 计算每个样本的SNP贡献
        std::vector<double> expected_snp_contributions(2, 0.0);
        for (size_t snp_idx = 0; snp_idx < 2; ++snp_idx)
        {
            const double p = frequencies[snp_idx];
            const double q = 1.0 - p;
            const double scale_add = std::sqrt(std::max(2.0 * p * q, eps));
            const double mean_add = 2.0 * p;
            const double effect = effects[snp_idx];

            for (size_t sample_idx = 0; sample_idx < 2; ++sample_idx)
            {
                const double x = geno(sample_idx, snp_idx);
                const double std_add = (x - mean_add) / scale_add;
                expected_snp_contributions[sample_idx] += std_add * effect;
            }
        }

        // 总预测值 = 截距 + SNP贡献
        std::vector<double> expected_predictions
            = {intercept + expected_snp_contributions[0],
               intercept + expected_snp_contributions[1]};

        using Catch::Matchers::WithinAbs;
        for (size_t i = 0; i < 2; ++i)
        {
            REQUIRE_THAT(
                predictions[i], WithinAbs(expected_predictions[i], eps));
            REQUIRE_THAT(
                snp_contributions[i],
                WithinAbs(expected_snp_contributions[i], eps));
            REQUIRE_THAT(covar_predictions(i, 0), WithinAbs(intercept, eps));
        }
    }

    SECTION("With quantitative covariate Age")
    {
        // 基因型与上一场景相同
        Eigen::MatrixXd genotypes(2, 2);
        genotypes << 0.0, 2.0, 1.0, 1.0;

        auto [bed_prefix, _] = bed_fixture.create_deterministic_bed_files(
            genotypes,
            {"sample1", "sample2"},
            {"rs1", "rs2"},
            {"1", "1"},
            {{'A', 'C'}, {'T', 'G'}});

        // 读取FAM文件获取FID/IID（用于协变量文件）
        auto fam_path = bed_prefix;
        fam_path.replace_extension(".fam");
        auto [fids, iids] = read_fam(fam_path);

        // SNP效应文件（同前）
        std::string snp_header = "Chrom\tPosition\tID\tA1\tA2\tA1Freq\tAdd";
        std::vector<std::string> snp_rows
            = {"1\t1000\trs1\tA\tC\t0.300000\t0.100000",
               "1\t2000\trs2\tT\tG\t0.400000\t-0.050000"};
        auto snp_path = file_fixture.create_text_file(
            create_snp_effect_content(snp_header, snp_rows), ".snp.eff");

        // 定量协变量文件：Age
        std::string qcovar_content = "FID\tIID\tAge\n";
        for (size_t i = 0; i < fids.size(); ++i)
        {
            qcovar_content += std::format(
                "{}\t{}\t{}\n", fids[i], iids[i], 25 + i * 5);  // 25, 30
        }
        auto qcovar_path
            = file_fixture.create_text_file(qcovar_content, ".qcovar");

        // 协变量效应文件：截距 1.0 + Age系数 0.2
        auto covar_path = file_fixture.create_text_file(
            "term\tmean\tstddev\tpercentile_5\tpercentile_95\tess\trhat\n"
            "Intercept\t1.0\t0.1\t0.8\t1.2\t1000\t1.0\n"
            "Age\t0.2\t0.05\t0.1\t0.3\t800\t1.01\n",
            ".param");

        auto output_path
            = file_fixture.generate_random_file_path(".predictions");

        PredictEngine::Config config;
        config.bed_path = bed_prefix;
        config.snp_effect_path = snp_path;
        config.covar_effect_path = covar_path;
        config.qcovar_path = qcovar_path;
        config.dcovar_path = "";
        config.output_path = output_path;
        config.iid_only = false;

        PredictEngine engine(config);
        REQUIRE_NOTHROW(engine.run());

        const auto& predictions = engine.predictions();
        const auto& snp_contributions = engine.snp_predictions();
        const auto& covar_predictions = engine.covar_predictions();

        REQUIRE(predictions.size() == 2);
        REQUIRE(snp_contributions.size() == 2);
        REQUIRE(covar_predictions.rows() == 2);
        REQUIRE(covar_predictions.cols() == 2);  // Intercept + Age

        // 精确计算期望值（使用与 PredictEngine 相同的公式）
        const double intercept = 1.0;
        const double age_coef = 0.2;
        const std::vector<double> ages = {25.0, 30.0};
        const double eps = 1e-10;

        // 精确计算SNP贡献（与上一场景相同）
        const std::vector<double> frequencies = {0.3, 0.4};
        const std::vector<double> effects = {0.1, -0.05};
        std::vector<double> expected_snp_contributions(2, 0.0);
        for (size_t snp_idx = 0; snp_idx < 2; ++snp_idx)
        {
            const double p = frequencies[snp_idx];
            const double q = 1.0 - p;
            const double scale_add = std::sqrt(std::max(2.0 * p * q, eps));
            const double mean_add = 2.0 * p;
            const double effect = effects[snp_idx];

            for (size_t sample_idx = 0; sample_idx < 2; ++sample_idx)
            {
                const double x = genotypes(sample_idx, snp_idx);
                const double std_add = (x - mean_add) / scale_add;
                expected_snp_contributions[sample_idx] += std_add * effect;
            }
        }

        using Catch::Matchers::WithinAbs;
        for (size_t i = 0; i < 2; ++i)
        {
            // 协变量贡献：intercept + age_coef * ages[i]
            double expected_covar_intercept = intercept;
            double expected_covar_age = age_coef * ages[i];
            double expected_total = expected_covar_intercept
                                    + expected_covar_age
                                    + expected_snp_contributions[i];

            REQUIRE_THAT(predictions[i], WithinAbs(expected_total, eps));
            REQUIRE_THAT(
                snp_contributions[i],
                WithinAbs(expected_snp_contributions[i], eps));
            REQUIRE_THAT(
                covar_predictions(i, 0),
                WithinAbs(expected_covar_intercept, eps));
            REQUIRE_THAT(
                covar_predictions(i, 1), WithinAbs(expected_covar_age, eps));
        }
    }
}
