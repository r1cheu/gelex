/*
 * Copyright 2026 RuLei Chen
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 * http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

#include <fmt/format.h>
#include <cstddef>
#include <filesystem>
#include <fstream>
#include <sstream>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

#include <Eigen/Core>
#include <catch2/catch_test_macros.hpp>

#include "file_fixture.h"
#include "gelex/data/genotype/process_method.h"
#include "gelex/data/genotype/processor.h"
#include "gelex/engine/predict.h"
#include "gelex/io/locistats/writer.h"
#include "gelex/types/genetic_effect_type.h"

#include "bed_fixture.h"

using gelex::PredictEngine;
using gelex::test::BedFixture;

namespace
{

// ---- file helpers -----------------------------------------------------------

auto write_tsv_row(const std::vector<std::string>& cols) -> std::string
{
    std::string line;
    for (size_t i = 0; i < cols.size(); ++i)
    {
        line += cols[i];
        if (i + 1 < cols.size())
        {
            line += '\t';
        }
    }
    line += '\n';
    return line;
}

auto make_tsv(
    const std::vector<std::string>& headers,
    const std::vector<std::vector<std::string>>& rows) -> std::string
{
    std::string out = write_tsv_row(headers);
    for (const auto& row : rows)
    {
        out += write_tsv_row(row);
    }
    return out;
}

auto read_fam(const std::filesystem::path& fam_path)
    -> std::pair<std::vector<std::string>, std::vector<std::string>>
{
    std::vector<std::string> fids;
    std::vector<std::string> iids;
    std::ifstream ifs(fam_path);
    std::string line;
    while (std::getline(ifs, line))
    {
        std::istringstream iss(line);
        std::string fid;
        std::string iid;
        iss >> fid >> iid;
        fids.push_back(std::move(fid));
        iids.push_back(std::move(iid));
    }
    return {std::move(fids), std::move(iids)};
}

// sbin 写入并返回标准化后基因型矩阵（add + dom）
auto create_sbin(
    const std::filesystem::path& sbin_path,
    const Eigen::MatrixXd& genotypes) -> void
{
    using gelex::EffectType;
    using gelex::GeneticMode;
    using gelex::GenotypeProcessMethod;
    using gelex::LociStatsWriter;
    using gelex::genotype::StandardizeHWE;

    const auto n_snps = genotypes.cols();

    Eigen::MatrixXd add_geno = genotypes;
    Eigen::VectorXd add_mean(n_snps);
    Eigen::VectorXd add_stddev(n_snps);
    for (Eigen::Index j = 0; j < n_snps; ++j)
    {
        auto col = add_geno.col(j);
        auto stats = StandardizeHWE<GeneticMode::A>::process(col);
        add_mean(j) = stats.mean;
        add_stddev(j) = stats.stddev;
    }

    Eigen::MatrixXd dom_geno = genotypes;
    Eigen::VectorXd dom_mean(n_snps);
    Eigen::VectorXd dom_stddev(n_snps);
    for (Eigen::Index j = 0; j < n_snps; ++j)
    {
        auto col = dom_geno.col(j);
        auto stats = StandardizeHWE<GeneticMode::D>::process(col);
        dom_mean(j) = stats.mean;
        dom_stddev(j) = stats.stddev;
    }

    LociStatsWriter writer(sbin_path.string());
    writer.write(
        EffectType::add(),
        GenotypeProcessMethod::StandardizeHWE().to_byte(),
        add_mean,
        &add_stddev);
    writer.write(
        EffectType::dom(),
        GenotypeProcessMethod::StandardizeHWE().to_byte(),
        dom_mean,
        &dom_stddev);
}

auto create_snp_effects_file(
    gelex::test::FileFixture& ff,
    std::string_view gfile_prefix,
    const std::vector<std::vector<std::string>>& snp_rows) -> void
{
    std::vector<std::string> headers
        = {"CHR", "SNP", "BP", "A1", "A2", "A1FREQ", "BETA_A", "BETA_D"};
    (void)ff.create_named_text_file(
        std::string(gfile_prefix) + ".snp.eff", make_tsv(headers, snp_rows));
}

auto create_param_file(
    gelex::test::FileFixture& ff,
    std::string_view gfile_prefix,
    double intercept,
    const std::vector<std::pair<std::string, double>>& qcovar_coefs,
    const std::vector<std::pair<std::string, double>>& dcovar_coefs) -> void
{
    const std::vector<std::string> param_headers
        = {"term",
           "mean",
           "stddev",
           "percentile_5",
           "percentile_95",
           "ess",
           "rhat"};

    std::vector<std::vector<std::string>> rows;
    rows.push_back(
        {"Intercept",
         std::to_string(intercept),
         "0.1",
         "0.8",
         "1.2",
         "1000",
         "1.0"});
    for (const auto& [name, coef] : qcovar_coefs)
    {
        rows.push_back(
            {name, std::to_string(coef), "0.05", "0.1", "0.3", "800", "1.01"});
    }
    for (const auto& [name, coef] : dcovar_coefs)
    {
        rows.push_back(
            {name,
             std::to_string(coef),
             "0.02",
             "-0.34",
             "-0.26",
             "900",
             "1.02"});
    }
    (void)ff.create_named_text_file(
        std::string(gfile_prefix) + ".param", make_tsv(param_headers, rows));
}

auto create_qcovar_file(
    gelex::test::FileFixture& ff,
    const std::vector<std::string>& fids,
    const std::vector<std::string>& iids,
    const std::vector<std::pair<std::string, std::vector<double>>>& covars)
    -> std::filesystem::path
{
    std::string content = "FID\tIID";
    for (const auto& [name, _] : covars)
    {
        content += '\t';
        content += name;
    }
    content += '\n';
    for (size_t i = 0; i < iids.size(); ++i)
    {
        content += fids[i] + '\t' + iids[i];
        for (const auto& [_, vals] : covars)
        {
            content += '\t';
            content += fmt::format("{}", vals[i]);
        }
        content += '\n';
    }
    return ff.create_text_file(content, ".qcovar");
}

auto create_dcovar_file(
    gelex::test::FileFixture& ff,
    const std::vector<std::string>& fids,
    const std::vector<std::string>& iids,
    const std::vector<std::pair<std::string, std::vector<std::string>>>& covars)
    -> std::filesystem::path
{
    std::string content = "FID\tIID";
    for (const auto& [name, _] : covars)
    {
        content += '\t';
        content += name;
    }
    content += '\n';
    for (size_t i = 0; i < iids.size(); ++i)
    {
        content += fids[i] + '\t' + iids[i];
        for (const auto& [_, vals] : covars)
        {
            content += '\t';
            content += vals[i];
        }
        content += '\n';
    }
    return ff.create_text_file(content, ".dcovar");
}

// ---- output reader ----------------------------------------------------------

struct PredictionOutput
{
    std::vector<std::string> headers;
    size_t row_count{};
    bool has_dominant{false};
};

auto read_prediction_output(const std::filesystem::path& path)
    -> PredictionOutput
{
    PredictionOutput out;
    std::ifstream ifs(path);
    std::string line;

    std::getline(ifs, line);
    {
        std::istringstream iss(line);
        std::string tok;
        while (std::getline(iss, tok, '\t'))
        {
            out.headers.push_back(tok);
        }
    }
    out.has_dominant = !out.headers.empty() && out.headers.back() == "dominant";

    while (std::getline(ifs, line))
    {
        ++out.row_count;
    }
    return out;
}

}  // namespace

TEST_CASE(
    "PredictEngine smoke test - full model",
    "[predict][predict_engine][smoke]")
{
    // 3 samples, 2 SNPs, add + dom + qcovar(Age) + dcovar(Sex)
    BedFixture bed;

    Eigen::MatrixXd genotypes(3, 2);
    genotypes << 0.0, 2.0, 1.0, 1.0, 2.0, 0.0;

    const std::vector<std::string> iids = {"s1", "s2", "s3"};
    const std::vector<std::string> snp_ids = {"rs1", "rs2"};
    const std::vector<std::pair<char, char>> alleles = {{'A', 'C'}, {'T', 'G'}};

    auto [bed_prefix, _] = bed.create_deterministic_bed_files(
        genotypes,
        iids,
        snp_ids,
        std::vector<std::string>(snp_ids.size(), "1"),
        alleles);

    auto [fids, loaded_iids] = read_fam(bed_prefix.string() + ".fam");

    auto& ff = bed.get_file_fixture();

    const std::vector<std::vector<std::string>> snp_rows
        = {{"1", "rs1", "1000", "A", "C", "0.30", "0.10", "0.02"},
           {"1", "rs2", "2000", "T", "G", "0.40", "-0.05", "0.01"}};

    auto gfile_prefix = (ff.get_test_dir() / "gfile").string();
    create_snp_effects_file(ff, gfile_prefix, snp_rows);
    create_param_file(
        ff, gfile_prefix, 1.0, {{"Age", 0.2}}, {{"Sex\x1FM", -0.3}});
    create_sbin(gfile_prefix + ".sbin", genotypes);

    auto qcovar_path = create_qcovar_file(
        ff, fids, loaded_iids, {{"Age", {25.0, 30.0, 35.0}}});
    auto dcovar_path
        = create_dcovar_file(ff, fids, loaded_iids, {{"Sex", {"M", "F", "M"}}});

    auto output_path = ff.get_test_dir() / "test.predictions";

    PredictEngine::Config config{
        .bfile_prefix = bed_prefix.string(),
        .gfile_prefix = gfile_prefix,
        .qcovar_path = qcovar_path,
        .dcovar_path = dcovar_path,
        .output_path = output_path};

    PredictEngine engine(config);
    REQUIRE_NOTHROW(engine.run());

    REQUIRE(std::filesystem::exists(output_path));

    auto out = read_prediction_output(output_path);
    REQUIRE(out.row_count == 3);
    REQUIRE(out.has_dominant);
}
