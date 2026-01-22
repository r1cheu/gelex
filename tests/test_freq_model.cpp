#include <filesystem>
#include <format>
#include <string>
#include <vector>

#include <Eigen/Core>
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include "../src/data/grm_bin_writer.h"
#include "../src/data/grm_id_writer.h"
#include "bed_fixture.h"
#include "file_fixture.h"
#include "gelex/data/data_pipe.h"
#include "gelex/model/freq/model.h"

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
    explicit GrmFileFixture(FileFixture& files, const std::string& prefix = "")
        : files_(files),
          prefix_(
              prefix.empty() ? files.generate_random_file_path("")
                             : std::filesystem::path(prefix))
    {
    }

    auto create(
        const Eigen::MatrixXd& matrix,
        const std::vector<std::string>& ids,
        double denominator = 1.0) -> void
    {
        auto bin_path = fs::path(prefix_.string() + ".grm.bin");
        {
            detail::GrmBinWriter writer(bin_path);
            writer.write(matrix, denominator);
        }

        auto id_path = fs::path(prefix_.string() + ".grm.id");
        {
            detail::GrmIdWriter writer(id_path);
            writer.write(ids);
        }
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

// Helper to create sample IDs in FID_IID format
auto make_sample_ids(Eigen::Index n, const std::string& prefix = "fam")
    -> std::vector<std::string>
{
    std::vector<std::string> ids;
    ids.reserve(static_cast<size_t>(n));
    for (Eigen::Index i = 0; i < n; ++i)
    {
        ids.push_back(std::format("{}{}_{}", prefix, i + 1, i + 1));
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
        // split "FID_IID" -> "FID", "IID"
        auto underscore_pos = ids[i].find('_');
        std::string fid = ids[i].substr(0, underscore_pos);
        std::string iid = ids[i].substr(underscore_pos + 1);
        content += std::format("{}\t{}\t{}\n", fid, iid, values(i));
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
        auto underscore_pos = ids[i].find('_');
        std::string fid = ids[i].substr(0, underscore_pos);
        std::string iid = ids[i].substr(underscore_pos + 1);
        content += std::format("{}\t{}", fid, iid);
        for (Eigen::Index j = 0; j < values.cols(); ++j)
        {
            content += std::format("\t{}", values(i, j));
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
        auto underscore_pos = ids[i].find('_');
        std::string fid = ids[i].substr(0, underscore_pos);
        std::string iid = ids[i].substr(underscore_pos + 1);
        content += std::format("{}\t{}", fid, iid);
        for (size_t j = 0; j < values[i].size(); ++j)
        {
            content += "\t" + values[i][j];
        }
        content += "\n";
    }
    return content;
}

}  // namespace

// ============================================================================
// FreqModel construction tests via DataPipe
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
    auto sample_ids = make_sample_ids(num_samples, "fam");

    // BedFixture creates samples as "fam{i%5+1}_sample{i+1}"
    std::vector<std::string> bed_sample_ids;
    bed_sample_ids.reserve(num_samples);
    for (Eigen::Index i = 0; i < num_samples; ++i)
    {
        bed_sample_ids.push_back(
            std::format("fam{}_sample{}", (i % 5) + 1, i + 1));
    }

    Eigen::VectorXd pheno_values
        = Eigen::VectorXd::LinSpaced(num_samples, 1, 10);
    auto pheno_content = make_phenotype_content(bed_sample_ids, pheno_values);
    auto pheno_path = files.create_text_file(pheno_content, ".phen");

    DataPipe::Config config{
        .phenotype_path = pheno_path,
        .phenotype_column = 2,  // 0-indexed: FID=0, IID=1, Phenotype=2
        .bed_path = bed_prefix,
    };

    DataPipe pipe(config);
    pipe.load_phenotypes();
    pipe.intersect_samples();
    pipe.finalize();

    FreqModel model(pipe);

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
        sample_ids.push_back(std::format("fam{}_sample{}", (i % 5) + 1, i + 1));
    }

    // create phenotype file
    Eigen::VectorXd pheno_values = Eigen::VectorXd::Random(num_samples);
    auto pheno_content = make_phenotype_content(sample_ids, pheno_values);
    auto pheno_path = files.create_text_file(pheno_content, ".phen");

    // create GRM files
    GrmFileFixture grm_fixture(files, "test.add");
    auto grm_matrix = make_symmetric_grm(num_samples);
    grm_fixture.create(grm_matrix, sample_ids);

    DataPipe::Config config{
        .phenotype_path = pheno_path,
        .phenotype_column = 2,  // 0-indexed
        .bed_path = bed_prefix,
        .grm_paths = {grm_fixture.prefix()}};

    DataPipe pipe(config);
    pipe.load_phenotypes();
    pipe.load_grms();
    pipe.intersect_samples();
    pipe.finalize();

    FreqModel model(pipe);

    SECTION("Verify one genetic effect")
    {
        REQUIRE(model.genetic().size() == 1);
    }

    SECTION("Verify genetic effect name is 'Additive'")
    {
        REQUIRE(model.genetic()[0].type == gelex::freq::GrmType::A);
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
        sample_ids.push_back(std::format("fam{}_sample{}", (i % 5) + 1, i + 1));
    }

    Eigen::VectorXd pheno_values = Eigen::VectorXd::Random(num_samples);
    auto pheno_content = make_phenotype_content(sample_ids, pheno_values);
    auto pheno_path = files.create_text_file(pheno_content, ".phen");

    GrmFileFixture grm_fixture(files, "test.dom");
    auto grm_matrix = make_symmetric_grm(num_samples);
    grm_fixture.create(grm_matrix, sample_ids);

    DataPipe::Config config{
        .phenotype_path = pheno_path,
        .phenotype_column = 2,  // 0-indexed
        .bed_path = bed_prefix,
        .grm_paths = {grm_fixture.prefix()}};

    DataPipe pipe(config);
    pipe.load_phenotypes();
    pipe.load_grms();
    pipe.intersect_samples();
    pipe.finalize();

    FreqModel model(pipe);

    SECTION("Verify one genetic effect")
    {
        REQUIRE(model.genetic().size() == 1);
    }

    SECTION("Verify genetic effect name is 'Dominance'")
    {
        REQUIRE(model.genetic()[0].type == gelex::freq::GrmType::D);
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
        sample_ids.push_back(std::format("fam{}_sample{}", (i % 5) + 1, i + 1));
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

    DataPipe::Config config{
        .phenotype_path = pheno_path,
        .phenotype_column = 2,  // 0-indexed
        .bed_path = bed_prefix,
        .grm_paths = {add_grm_fixture.prefix(), dom_grm_fixture.prefix()}};

    DataPipe pipe(config);
    pipe.load_phenotypes();
    pipe.load_grms();
    pipe.intersect_samples();
    pipe.finalize();

    FreqModel model(pipe);

    SECTION("Verify two genetic effects")
    {
        REQUIRE(model.genetic().size() == 2);
    }

    SECTION("Verify genetic effect names")
    {
        REQUIRE(model.genetic()[0].type == gelex::freq::GrmType::A);
        REQUIRE(model.genetic()[1].type == gelex::freq::GrmType::D);
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
        bed_sample_ids.push_back(
            std::format("fam{}_sample{}", (i % 5) + 1, i + 1));
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

    DataPipe::Config config{
        .phenotype_path = pheno_path,
        .phenotype_column = 2,  // 0-indexed
        .bed_path = bed_prefix,
        .grm_paths = {grm_fixture.prefix()}};

    DataPipe pipe(config);
    pipe.load_phenotypes();
    pipe.load_grms();
    auto stats = pipe.intersect_samples();
    pipe.finalize();

    FreqModel model(pipe);

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
        sample_ids.push_back(std::format("fam{}_sample{}", (i % 5) + 1, i + 1));
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

    DataPipe::Config config{
        .phenotype_path = pheno_path,
        .phenotype_column = 2,  // 0-indexed
        .bed_path = bed_prefix,
        .qcovar_path = qcovar_path,
        .grm_paths = {grm_fixture.prefix()}};

    DataPipe pipe(config);
    pipe.load_phenotypes();
    pipe.load_covariates();
    pipe.load_grms();
    pipe.intersect_samples();
    pipe.finalize();

    FreqModel model(pipe);

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
        sample_ids.push_back(std::format("fam{}_sample{}", (i % 5) + 1, i + 1));
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

    DataPipe::Config config{
        .phenotype_path = pheno_path,
        .phenotype_column = 2,  // 0-indexed
        .bed_path = bed_prefix,
        .dcovar_path = dcovar_path,
        .grm_paths = {grm_fixture.prefix()}};

    DataPipe pipe(config);
    pipe.load_phenotypes();
    pipe.load_covariates();
    pipe.load_grms();
    pipe.intersect_samples();
    pipe.finalize();

    FreqModel model(pipe);

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
        sample_ids.push_back(std::format("fam{}_sample{}", (i % 5) + 1, i + 1));
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

    DataPipe::Config config{
        .phenotype_path = pheno_path,
        .phenotype_column = 2,  // 0-indexed
        .bed_path = bed_prefix,
        .grm_paths = {grm_fixture.prefix()}};

    DataPipe pipe(config);
    pipe.load_phenotypes();
    pipe.load_grms();
    pipe.intersect_samples();
    pipe.finalize();

    FreqModel model(pipe);

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
