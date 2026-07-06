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
#include <optional>
#include <sstream>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

#include <Eigen/Core>
#include <catch2/catch_test_macros.hpp>

#include "file_fixture.h"
#include "gelex/data/bed.h"
#include "gelex/data/dataframe/dataframe.h"
#include "gelex/data/genotype_method.h"
#include "gelex/data/locus_encoding.h"
#include "gelex/data/reader.h"
#include "gelex/exception.h"
#include "gelex/io/binary_reader.h"
#include "gelex/io/predict_reader.h"
#include "gelex/io/predict_writer.h"
#include "gelex/io/snpstats.h"
#include "gelex/predict/compute.h"
#include "gelex/predict/snp_alignment.h"
#include "gelex/predict/standardize.h"
#include "gelex/predict/types.h"
#include "gelex/types/genetic_mode.h"

#include "bed_fixture.h"

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

// snpstats 写入并返回标准化后基因型矩阵（add + dom）
auto create_snpstats(
    const std::filesystem::path& snpstats_path,
    const Eigen::MatrixXd& genotypes) -> void
{
    using gelex::BinaryWriter;
    using gelex::GeneticMode;
    using gelex::GenotypeMethod;
    using gelex::SnpStats;
    using gelex::write_snp_stats;

    const auto n_snps = genotypes.cols();

    const auto fill = [&genotypes, n_snps](GeneticMode mode) -> SnpStats
    {
        const gelex::EncodingSpec spec{gelex::encoding_spec_from_method(
            mode, GenotypeMethod::StandardizeHWE)};
        const gelex::LociEncoding encoding{
            gelex::detail::make_loci_encoding<double>(genotypes, spec)};

        SnpStats stats;
        stats.method = GenotypeMethod::StandardizeHWE;
        stats.mean = Eigen::VectorXd(n_snps);
        stats.var = Eigen::VectorXd(n_snps);
        stats.A1freq = Eigen::VectorXd(n_snps);
        for (const gelex::LocusEncoding& locus : encoding.loci)
        {
            stats.mean(locus.marker_index) = locus.mean;
            stats.var(locus.marker_index) = locus.var;
            stats.A1freq(locus.marker_index)
                = locus.stats.has_nonmissing() ? locus.stats.A1freq() : 0.0;
        }
        return stats;
    };

    BinaryWriter writer(snpstats_path.string());
    write_snp_stats(writer, GeneticMode::A, fill(GeneticMode::A));
    write_snp_stats(writer, GeneticMode::D, fill(GeneticMode::D));
}

auto create_snp_effects_file(
    gelex::test::FileFixture& ff,
    std::string_view gfile_prefix,
    const std::vector<std::vector<std::string>>& snp_rows) -> void
{
    std::vector<std::string> headers
        = {"CHR", "SNP", "BP", "A1", "A2", "A1FREQ", "BETA_A", "BETA_D"};
    (void)ff.create_named_text_file(
        std::string(gfile_prefix) + ".snpeff", make_tsv(headers, snp_rows));
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

auto load_snpstats(const std::filesystem::path& path) -> gelex::SnpStatsData
{
    gelex::BinaryReader reader(path.string());
    gelex::SnpStatsData data;
    data.add = gelex::read_snp_stats(reader, gelex::GeneticMode::A);
    if (gelex::has_snp_stats(reader, gelex::GeneticMode::D))
    {
        data.dom = gelex::read_snp_stats(reader, gelex::GeneticMode::D);
        data.has_dom = true;
    }
    return data;
}

auto run_predict_dataflow(
    const std::string& bfile_prefix,
    const std::string& gfile_prefix,
    const std::optional<std::filesystem::path>& qcovar_path,
    const std::optional<std::filesystem::path>& dcovar_path,
    const std::filesystem::path& output_path) -> void
{
    auto snp_effects = gelex::read_snp_effects(gfile_prefix + ".snpeff");
    auto snpstats = load_snpstats(gfile_prefix + ".snpstats");

    bool enable_dom{};
    if (snpstats.has_dom)
    {
        if (!snp_effects.contains("BETA_D"))
        {
            throw gelex::GelexException(
                ".snpstats file contains dominance effects, but SNP effects "
                "file "
                "does not have 'BETA_D' column.");
        }
        enable_dom = true;
    }

    Eigen::VectorXd add_effects = snp_effects["BETA_A"].to_map<double>();
    auto dom_effects = enable_dom ? std::make_optional<Eigen::VectorXd>(
                                        snp_effects["BETA_D"].to_mat<double>())
                                  : std::nullopt;
    auto coefficients = gelex::read_coefficients(gfile_prefix + ".param");

    auto fam_df = gelex::read_fam(bfile_prefix + ".fam");
    auto bim_df = gelex::read_bim(bfile_prefix + ".bim");
    auto covariates = gelex::read_covariates(
        qcovar_path, dcovar_path, coefficients, fam_df);

    auto alignment = gelex::build_snp_alignment(snp_effects, bim_df);
    if (alignment.num_missing > 0 || alignment.num_mismatched > 0)
    {
        throw gelex::GelexException(
            fmt::format(
                "{}.snpeff does not match {}.bim: {} missing SNPs, {} "
                "allele mismatches",
                gfile_prefix,
                bfile_prefix,
                alignment.num_missing,
                alignment.num_mismatched));
    }

    auto bed = gelex::open_bed(bfile_prefix, fam_df.index());
    auto genotype = bed.read_snps<double>(alignment.column_map);

    gelex::GenotypeData geno;
    if (snpstats.has_dom)
    {
        geno.dom = genotype;
    }
    geno.add = std::move(genotype);

    gelex::standardize_genotypes(geno, snpstats);

    gelex::SnpEffects effects{
        .add = std::move(add_effects), .dom = std::move(dom_effects)};
    auto gebv = gelex::compute_gebv(geno, effects);
    auto covar = gelex::compute_covariate_effects(covariates, coefficients);

    auto sample_keys = fam_df.index().keys();
    std::vector<std::string> sample_ids(sample_keys.begin(), sample_keys.end());

    gelex::PredictResult result{
        .sample_ids = std::move(sample_ids),
        .predictions = gebv.total + covar.total,
        .snp_predictions = std::move(gebv.total),
        .add_predictions = std::move(gebv.add_predictions),
        .dom_predictions = std::move(gebv.dom_predictions),
        .covar_predictions = std::move(covar.per_covariate),
        .covar_names = std::move(covar.covar_names)};

    gelex::PredictWriter writer(output_path);
    writer.write(result);
}

}  // namespace

TEST_CASE(
    "Predict command dataflow writes full model predictions",
    "[predict][dataflow][smoke]")
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
    create_snpstats(gfile_prefix + ".snpstats", genotypes);

    auto qcovar_path = create_qcovar_file(
        ff, fids, loaded_iids, {{"Age", {25.0, 30.0, 35.0}}});
    auto dcovar_path
        = create_dcovar_file(ff, fids, loaded_iids, {{"Sex", {"M", "F", "M"}}});

    auto output_path = ff.get_test_dir() / "test.predictions";

    REQUIRE_NOTHROW(run_predict_dataflow(
        bed_prefix.string(),
        gfile_prefix,
        qcovar_path,
        dcovar_path,
        output_path));

    REQUIRE(std::filesystem::exists(output_path));

    auto out = read_prediction_output(output_path);
    REQUIRE(out.row_count == 3);
    REQUIRE(out.has_dominant);
}

TEST_CASE(
    "Predict command dataflow rejects missing SNPs",
    "[predict][dataflow]")
{
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

    auto& ff = bed.get_file_fixture();

    const std::vector<std::vector<std::string>> snp_rows
        = {{"1", "rs1", "1000", "A", "C", "0.30", "0.10", "0.02"},
           {"1", "rs3", "3000", "A", "G", "0.40", "-0.05", "0.01"}};

    auto gfile_prefix = (ff.get_test_dir() / "missing_snp").string();
    create_snp_effects_file(ff, gfile_prefix, snp_rows);
    create_param_file(ff, gfile_prefix, 1.0, {}, {});
    create_snpstats(gfile_prefix + ".snpstats", genotypes);

    auto output_path = ff.get_test_dir() / "missing.predictions";

    REQUIRE_THROWS_AS(
        run_predict_dataflow(
            bed_prefix.string(),
            gfile_prefix,
            std::nullopt,
            std::nullopt,
            output_path),
        gelex::GelexException);
}
