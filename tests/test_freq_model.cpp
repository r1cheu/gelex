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
#include <string>
#include <utility>
#include <vector>

#include <Eigen/Core>
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include "bed_fixture.h"
#include "file_fixture.h"
#include "gelex/data/dataframe/index.h"
#include "gelex/data/pipe/grm.h"
#include "gelex/data/pipe/pheno.h"
#include "gelex/data/reader.h"
#include "gelex/data/sample_id.h"
#include "gelex/io/grm/writer.h"
#include "gelex/model/freq/model.h"
#include "gelex/types/genetic_effect_type.h"
#include "sample_id_fixture.h"

namespace fs = std::filesystem;

using namespace gelex;        // NOLINT
using namespace gelex::test;  // NOLINT
using Catch::Matchers::WithinAbs;

namespace
{

// Helper class to create GRM test files
class GrmFileFixture
{
   public:
    explicit GrmFileFixture(FileFixture& files, const std::string& suffix = "")
        : files_(files),
          prefix_(
              suffix.empty() ? files.generate_random_file_path("")
                             : files.get_test_dir() / suffix)
    {
    }

    auto create(
        const Eigen::MatrixXd& matrix,
        const std::vector<std::string>& ids) -> void
    {
        auto bin_path = fs::path(prefix_.string() + ".bin");
        {
            gelex::GrmBinWriter writer(bin_path);
            writer.write(matrix);
        }

        gelex::write_grm_ids(prefix_.string() + ".id", ids);
    }

    [[nodiscard]] auto prefix() const -> const fs::path& { return prefix_; }

   private:
    FileFixture& files_;
    fs::path prefix_;
};

// Helper to create a symmetric GRM matrix
auto make_symmetric_grm(Eigen::Index n) -> Eigen::MatrixXd
{
    Eigen::MatrixXd matrix = Eigen::MatrixXd::Random(n, n);
    matrix = (matrix + matrix.transpose()) / 2.0;
    // ensure diagonal is positive (like a real GRM)
    matrix.diagonal().array() += 2.0;
    return matrix;
}

auto sid(std::string_view fid, std::string_view iid) -> std::string
{
    return gelex::make_sample_id(fid, iid);
}

// Helper to create sample IDs in canonical format
auto make_sample_ids(Eigen::Index n, const std::string& prefix = "fam")
    -> std::vector<std::string>
{
    std::vector<std::string> ids;
    ids.reserve(static_cast<size_t>(n));
    for (Eigen::Index i = 0; i < n; ++i)
    {
        ids.push_back(
            sid(fmt::format("{}{}", prefix, i + 1), std::to_string(i + 1)));
    }
    return ids;
}

// Helper to create phenotype file content
auto make_phenotype_content(
    const std::vector<std::string>& ids,
    const Eigen::VectorXd& values) -> std::string
{
    std::string content = "FID\tIID\tPhenotype\n";
    for (size_t i = 0; i < ids.size(); ++i)
    {
        auto [fid, iid] = gelex::split_sample_id(ids[i]);
        content += fmt::format("{}\t{}\t{}\n", fid, iid, values(i));
    }
    return content;
}

// Helper to create qcovar file content
auto make_qcovar_content(
    const std::vector<std::string>& ids,
    const Eigen::MatrixXd& values,
    const std::vector<std::string>& col_names) -> std::string
{
    std::string content = "FID\tIID";
    for (const auto& name : col_names)
    {
        content += "\t" + name;
    }
    content += "\n";

    for (size_t i = 0; i < ids.size(); ++i)
    {
        auto [fid, iid] = gelex::split_sample_id(ids[i]);
        content += fmt::format("{}\t{}", fid, iid);
        for (Eigen::Index j = 0; j < values.cols(); ++j)
        {
            content += fmt::format("\t{}", values(i, j));
        }
        content += "\n";
    }
    return content;
}

// Helper to create dcovar file content
auto make_dcovar_content(
    const std::vector<std::string>& ids,
    const std::vector<std::vector<std::string>>& values,
    const std::vector<std::string>& col_names) -> std::string
{
    std::string content = "FID\tIID";
    for (const auto& name : col_names)
    {
        content += "\t" + name;
    }
    content += "\n";

    for (size_t i = 0; i < ids.size(); ++i)
    {
        auto [fid, iid] = gelex::split_sample_id(ids[i]);
        content += fmt::format("{}\t{}", fid, iid);
        for (size_t j = 0; j < values[i].size(); ++j)
        {
            content += "\t" + values[i][j];
        }
        content += "\n";
    }
    return content;
}

// Helper to build PhenoPipe + GrmPipe and FreqModel
auto make_freq_model(
    PhenoPipe::Config pheno_config,
    const fs::path& bed_prefix,
    const std::vector<std::filesystem::path>& grm_paths = {}) -> FreqModel
{
    auto fam_index = read_fam(bed_prefix.string() + ".fam").index();

    PhenoPipe pheno(pheno_config);

    GrmPipe grm(grm_paths);

    std::vector<const dataframe::Index<std::string>*> all_indices{&fam_index};
    all_indices.append_range(pheno.sample_indices());
    all_indices.append_range(grm.sample_indices());
    auto common = dataframe::intersect<std::string>(all_indices);

    pheno.load(common);
    grm.load(common);

    return FreqModel(
        std::move(pheno).take_phenotype(),
        std::move(pheno).take_fixed_design(),
        std::move(grm).take_grms());
}

}  // namespace

// ============================================================================
// FreqModel construction tests via PhenoPipe + GrmPipe
// ============================================================================

TEST_CASE(
    "FreqModel - Construction with phenotype only (no GRM)",
    "[freq_model][integration]")
{
    BedFixture bed_fixture;
    const Eigen::Index num_samples = 10;
    const Eigen::Index num_snps = 5;

    auto [bed_prefix, _] = bed_fixture.create_bed_files(num_samples, num_snps);

    // create phenotype file with same samples
    auto& files = bed_fixture.get_file_fixture();

    // BedFixture creates samples as "fam{i%5+1}_sample{i+1}"
    std::vector<std::string> bed_sample_ids;
    bed_sample_ids.reserve(num_samples);
    for (Eigen::Index i = 0; i < num_samples; ++i)
    {
        bed_sample_ids.push_back(sid(
            fmt::format("fam{}", (i % 5) + 1), fmt::format("sample{}", i + 1)));
    }

    Eigen::VectorXd pheno_values
        = Eigen::VectorXd::LinSpaced(num_samples, 1, 10);
    auto pheno_content = make_phenotype_content(bed_sample_ids, pheno_values);
    auto pheno_path = files.create_text_file(pheno_content, ".phen");

    PhenoPipe::Config pheno_config{
        .phenotype_path = pheno_path,
        .phenotype_column = 2,
    };

    auto model = make_freq_model(pheno_config, bed_prefix);

    SECTION("Verify num_individuals")
    {
        REQUIRE(model.num_individuals() == num_samples);
    }

    SECTION("Verify phenotype")
    {
        REQUIRE(model.phenotype().size() == num_samples);
    }

    SECTION("Verify fixed effects (intercept only)")
    {
        // without covariates, should have intercept column only
        REQUIRE(model.fixed().X.rows() == num_samples);
        REQUIRE(model.fixed().X.cols() >= 1);
    }

    SECTION("Verify no genetic effects")
    {
        REQUIRE(model.genetic().empty());
    }

    SECTION("Verify no random effects")
    {
        REQUIRE(model.random().empty());
    }
}

TEST_CASE(
    "FreqModel - Construction with additive GRM",
    "[freq_model][integration]")
{
    BedFixture bed_fixture;
    const Eigen::Index num_samples = 8;
    const Eigen::Index num_snps = 3;

    auto [bed_prefix, _] = bed_fixture.create_bed_files(num_samples, num_snps);
    auto& files = bed_fixture.get_file_fixture();

    // create sample IDs matching BedFixture format
    std::vector<std::string> sample_ids;
    sample_ids.reserve(num_samples);
    for (Eigen::Index i = 0; i < num_samples; ++i)
    {
        sample_ids.push_back(sid(
            fmt::format("fam{}", (i % 5) + 1), fmt::format("sample{}", i + 1)));
    }

    // create phenotype file
    Eigen::VectorXd pheno_values = Eigen::VectorXd::Random(num_samples);
    auto pheno_content = make_phenotype_content(sample_ids, pheno_values);
    auto pheno_path = files.create_text_file(pheno_content, ".phen");

    // create GRM files
    GrmFileFixture grm_fixture(files, "test.add");
    auto grm_matrix = make_symmetric_grm(num_samples);
    grm_fixture.create(grm_matrix, sample_ids);

    PhenoPipe::Config pheno_config{
        .phenotype_path = pheno_path,
        .phenotype_column = 2,
    };

    auto model
        = make_freq_model(pheno_config, bed_prefix, {grm_fixture.prefix()});

    SECTION("Verify one genetic effect")
    {
        REQUIRE(model.genetic().size() == 1);
    }

    SECTION("Verify genetic effect name is 'Additive'")
    {
        REQUIRE(model.genetic()[0].type == gelex::GeneticMode::A);
    }

    SECTION("Verify GRM matrix dimensions")
    {
        REQUIRE(model.genetic()[0].K.rows() == num_samples);
        REQUIRE(model.genetic()[0].K.cols() == num_samples);
    }

    SECTION("Verify GRM matrix is symmetric")
    {
        const auto& K = model.genetic()[0].K;
        for (Eigen::Index i = 0; i < K.rows(); ++i)
        {
            for (Eigen::Index j = 0; j < i; ++j)
            {
                REQUIRE(K(i, j) == K(j, i));
            }
        }
    }
}

TEST_CASE(
    "FreqModel - Construction with dominance GRM",
    "[freq_model][integration]")
{
    BedFixture bed_fixture;
    const Eigen::Index num_samples = 6;
    const Eigen::Index num_snps = 2;

    auto [bed_prefix, _] = bed_fixture.create_bed_files(num_samples, num_snps);
    auto& files = bed_fixture.get_file_fixture();

    std::vector<std::string> sample_ids;
    sample_ids.reserve(num_samples);
    for (Eigen::Index i = 0; i < num_samples; ++i)
    {
        sample_ids.push_back(sid(
            fmt::format("fam{}", (i % 5) + 1), fmt::format("sample{}", i + 1)));
    }

    Eigen::VectorXd pheno_values = Eigen::VectorXd::Random(num_samples);
    auto pheno_content = make_phenotype_content(sample_ids, pheno_values);
    auto pheno_path = files.create_text_file(pheno_content, ".phen");

    GrmFileFixture grm_fixture(files, "test.dom");
    auto grm_matrix = make_symmetric_grm(num_samples);
    grm_fixture.create(grm_matrix, sample_ids);

    PhenoPipe::Config pheno_config{
        .phenotype_path = pheno_path,
        .phenotype_column = 2,
    };

    auto model
        = make_freq_model(pheno_config, bed_prefix, {grm_fixture.prefix()});

    SECTION("Verify one genetic effect")
    {
        REQUIRE(model.genetic().size() == 1);
    }

    SECTION("Verify genetic effect name is 'Dominance'")
    {
        REQUIRE(model.genetic()[0].type == gelex::GeneticMode::D);
    }

    SECTION("Verify GRM matrix dimensions")
    {
        REQUIRE(model.genetic()[0].K.rows() == num_samples);
        REQUIRE(model.genetic()[0].K.cols() == num_samples);
    }
}

TEST_CASE(
    "FreqModel - Construction with both additive and dominance GRM",
    "[freq_model][integration]")
{
    BedFixture bed_fixture;
    const Eigen::Index num_samples = 5;
    const Eigen::Index num_snps = 2;

    auto [bed_prefix, _] = bed_fixture.create_bed_files(num_samples, num_snps);
    auto& files = bed_fixture.get_file_fixture();

    std::vector<std::string> sample_ids;
    sample_ids.reserve(num_samples);
    for (Eigen::Index i = 0; i < num_samples; ++i)
    {
        sample_ids.push_back(sid(
            fmt::format("fam{}", (i % 5) + 1), fmt::format("sample{}", i + 1)));
    }

    Eigen::VectorXd pheno_values = Eigen::VectorXd::Random(num_samples);
    auto pheno_content = make_phenotype_content(sample_ids, pheno_values);
    auto pheno_path = files.create_text_file(pheno_content, ".phen");

    // create additive GRM
    GrmFileFixture add_grm_fixture(files, "test.add");
    auto add_grm_matrix = make_symmetric_grm(num_samples);
    add_grm_fixture.create(add_grm_matrix, sample_ids);

    // create dominance GRM
    GrmFileFixture dom_grm_fixture(files, "test.dom");
    auto dom_grm_matrix = make_symmetric_grm(num_samples);
    dom_grm_fixture.create(dom_grm_matrix, sample_ids);

    PhenoPipe::Config pheno_config{
        .phenotype_path = pheno_path,
        .phenotype_column = 2,
    };

    auto model = make_freq_model(
        pheno_config,
        bed_prefix,
        {add_grm_fixture.prefix(), dom_grm_fixture.prefix()});

    SECTION("Verify two genetic effects")
    {
        REQUIRE(model.genetic().size() == 2);
    }

    SECTION("Verify genetic effect names")
    {
        REQUIRE(model.genetic()[0].type == gelex::GeneticMode::A);
        REQUIRE(model.genetic()[1].type == gelex::GeneticMode::D);
    }

    SECTION("Verify both GRM matrices have correct dimensions")
    {
        REQUIRE(model.genetic()[0].K.rows() == num_samples);
        REQUIRE(model.genetic()[0].K.cols() == num_samples);
        REQUIRE(model.genetic()[1].K.rows() == num_samples);
        REQUIRE(model.genetic()[1].K.cols() == num_samples);
    }
}

TEST_CASE(
    "FreqModel - Sample intersection with GRM",
    "[freq_model][integration]")
{
    BedFixture bed_fixture;
    const Eigen::Index bed_samples = 10;
    const Eigen::Index grm_samples = 8;
    const Eigen::Index num_snps = 2;

    auto [bed_prefix, _] = bed_fixture.create_bed_files(bed_samples, num_snps);
    auto& files = bed_fixture.get_file_fixture();

    // BED file has samples: fam1_sample1, fam2_sample2, ..., fam5_sample5,
    //                       fam1_sample6, fam2_sample7, ...
    std::vector<std::string> bed_sample_ids;
    bed_sample_ids.reserve(bed_samples);
    for (Eigen::Index i = 0; i < bed_samples; ++i)
    {
        bed_sample_ids.push_back(sid(
            fmt::format("fam{}", (i % 5) + 1), fmt::format("sample{}", i + 1)));
    }

    // GRM has only first 8 samples (subset)
    std::vector<std::string> grm_sample_ids(
        bed_sample_ids.begin(), bed_sample_ids.begin() + grm_samples);

    // phenotype has all BED samples
    Eigen::VectorXd pheno_values
        = Eigen::VectorXd::LinSpaced(bed_samples, 1, 10);
    auto pheno_content = make_phenotype_content(bed_sample_ids, pheno_values);
    auto pheno_path = files.create_text_file(pheno_content, ".phen");

    // GRM with subset of samples
    GrmFileFixture grm_fixture(files);
    auto grm_matrix = make_symmetric_grm(grm_samples);
    grm_fixture.create(grm_matrix, grm_sample_ids);

    PhenoPipe::Config pheno_config{
        .phenotype_path = pheno_path,
        .phenotype_column = 2,
    };

    auto model
        = make_freq_model(pheno_config, bed_prefix, {grm_fixture.prefix()});

    SECTION("Verify sample count reflects intersection")
    {
        // common samples = min(bed, grm, pheno) = 8
        REQUIRE(model.num_individuals() == grm_samples);
    }

    SECTION("Verify phenotype is filtered")
    {
        REQUIRE(model.phenotype().size() == grm_samples);
    }

    SECTION("Verify GRM matrix is filtered")
    {
        REQUIRE(model.genetic()[0].K.rows() == grm_samples);
        REQUIRE(model.genetic()[0].K.cols() == grm_samples);
    }

    SECTION("Verify fixed effects matrix is filtered")
    {
        REQUIRE(model.fixed().X.rows() == grm_samples);
    }
}

TEST_CASE(
    "FreqModel - Construction with quantitative covariates",
    "[freq_model][integration]")
{
    BedFixture bed_fixture;
    const Eigen::Index num_samples = 6;
    const Eigen::Index num_snps = 2;

    auto [bed_prefix, _] = bed_fixture.create_bed_files(num_samples, num_snps);
    auto& files = bed_fixture.get_file_fixture();

    std::vector<std::string> sample_ids;
    sample_ids.reserve(num_samples);
    for (Eigen::Index i = 0; i < num_samples; ++i)
    {
        sample_ids.push_back(sid(
            fmt::format("fam{}", (i % 5) + 1), fmt::format("sample{}", i + 1)));
    }

    // create phenotype
    Eigen::VectorXd pheno_values = Eigen::VectorXd::Random(num_samples);
    auto pheno_content = make_phenotype_content(sample_ids, pheno_values);
    auto pheno_path = files.create_text_file(pheno_content, ".phen");

    // create qcovar with 2 columns
    Eigen::MatrixXd qcovar_values(num_samples, 2);
    qcovar_values.col(0)
        = Eigen::VectorXd::LinSpaced(num_samples, 20, 50);  // age
    qcovar_values.col(1)
        = Eigen::VectorXd::LinSpaced(num_samples, 160, 185);  // height
    auto qcovar_content
        = make_qcovar_content(sample_ids, qcovar_values, {"Age", "Height"});
    auto qcovar_path = files.create_text_file(qcovar_content, ".qcovar");

    // create GRM
    GrmFileFixture grm_fixture(files);
    auto grm_matrix = make_symmetric_grm(num_samples);
    grm_fixture.create(grm_matrix, sample_ids);

    PhenoPipe::Config pheno_config{
        .phenotype_path = pheno_path,
        .phenotype_column = 2,
        .quantitative_covariates_path = qcovar_path,
    };

    auto model
        = make_freq_model(pheno_config, bed_prefix, {grm_fixture.prefix()});

    SECTION("Verify fixed effects include covariates")
    {
        // expected columns: intercept + 2 qcovars = 3
        REQUIRE(model.fixed().X.rows() == num_samples);
        REQUIRE(model.fixed().X.cols() == 3);
    }

    SECTION("Verify fixed effect names")
    {
        REQUIRE(model.fixed().names.size() == 3);
        REQUIRE(model.fixed().names[0] == "Intercept");
        REQUIRE(model.fixed().names[1] == "Age");
        REQUIRE(model.fixed().names[2] == "Height");
    }
}

TEST_CASE(
    "FreqModel - Construction with discrete covariates",
    "[freq_model][integration]")
{
    BedFixture bed_fixture;
    const Eigen::Index num_samples = 6;
    const Eigen::Index num_snps = 2;

    auto [bed_prefix, _] = bed_fixture.create_bed_files(num_samples, num_snps);
    auto& files = bed_fixture.get_file_fixture();

    std::vector<std::string> sample_ids;
    sample_ids.reserve(num_samples);
    for (Eigen::Index i = 0; i < num_samples; ++i)
    {
        sample_ids.push_back(sid(
            fmt::format("fam{}", (i % 5) + 1), fmt::format("sample{}", i + 1)));
    }

    // create phenotype
    Eigen::VectorXd pheno_values = Eigen::VectorXd::Random(num_samples);
    auto pheno_content = make_phenotype_content(sample_ids, pheno_values);
    auto pheno_path = files.create_text_file(pheno_content, ".phen");

    // create dcovar with 1 column (3 levels: A, B, C)
    std::vector<std::vector<std::string>> dcovar_values
        = {{"A"}, {"B"}, {"C"}, {"A"}, {"B"}, {"C"}};
    auto dcovar_content
        = make_dcovar_content(sample_ids, dcovar_values, {"Group"});
    auto dcovar_path = files.create_text_file(dcovar_content, ".dcovar");

    // create GRM
    GrmFileFixture grm_fixture(files);
    auto grm_matrix = make_symmetric_grm(num_samples);
    grm_fixture.create(grm_matrix, sample_ids);

    PhenoPipe::Config pheno_config{
        .phenotype_path = pheno_path,
        .phenotype_column = 2,
        .discrete_covariates_path = dcovar_path,
    };

    auto model
        = make_freq_model(pheno_config, bed_prefix, {grm_fixture.prefix()});

    SECTION("Verify fixed effects include dummy coded discrete covariate")
    {
        // expected columns: intercept + (3 levels - 1 reference) = 3
        REQUIRE(model.fixed().X.rows() == num_samples);
        REQUIRE(model.fixed().X.cols() == 3);
    }
}

TEST_CASE(
    "FreqModel - GRM values preserved after filtering",
    "[freq_model][integration]")
{
    BedFixture bed_fixture;
    const Eigen::Index num_samples = 4;
    const Eigen::Index num_snps = 2;

    auto [bed_prefix, _] = bed_fixture.create_bed_files(num_samples, num_snps);
    auto& files = bed_fixture.get_file_fixture();

    std::vector<std::string> sample_ids;
    sample_ids.reserve(num_samples);
    for (Eigen::Index i = 0; i < num_samples; ++i)
    {
        sample_ids.push_back(sid(
            fmt::format("fam{}", (i % 5) + 1), fmt::format("sample{}", i + 1)));
    }

    Eigen::VectorXd pheno_values = Eigen::VectorXd::Random(num_samples);
    auto pheno_content = make_phenotype_content(sample_ids, pheno_values);
    auto pheno_path = files.create_text_file(pheno_content, ".phen");

    // create known GRM matrix
    Eigen::MatrixXd original_grm(num_samples, num_samples);
    // clang-format off
    original_grm << 1.0, 0.5, 0.3, 0.2,
                    0.5, 1.0, 0.4, 0.1,
                    0.3, 0.4, 1.0, 0.6,
                    0.2, 0.1, 0.6, 1.0;
    // clang-format on

    GrmFileFixture grm_fixture(files);
    grm_fixture.create(original_grm, sample_ids);

    PhenoPipe::Config pheno_config{
        .phenotype_path = pheno_path,
        .phenotype_column = 2,
    };

    auto model
        = make_freq_model(pheno_config, bed_prefix, {grm_fixture.prefix()});

    SECTION("Verify GRM values match original (accounting for float precision)")
    {
        const auto& K = model.genetic()[0].K;
        for (Eigen::Index i = 0; i < num_samples; ++i)
        {
            for (Eigen::Index j = 0; j < num_samples; ++j)
            {
                // GRM is stored as float32, so we compare with float precision
                auto expected = static_cast<double>(
                    static_cast<float>(original_grm(i, j)));
                REQUIRE(K(i, j) == expected);
            }
        }
    }
}
