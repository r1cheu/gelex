#include <fstream>
#include <map>
#include <sstream>
#include <string>
#include <string_view>
#include <vector>

#include <Eigen/Core>
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_exception.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>

#include "../src/predict/predict_pipe.h"
#include "bed_fixture.h"
#include "file_fixture.h"

namespace fs = std::filesystem;

using namespace gelex;  // NOLINT
using Catch::Matchers::EndsWith;
using gelex::test::are_matrices_equal;
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
}  // namespace

TEST_CASE(
    "PredictDataPipe - Construction with BED only",
    "[predict][predict_pipe]")
{
    BedFixture bed_fixture;

    SECTION("Happy path - only BED file")
    {
        const Eigen::Index num_samples = 5;
        const Eigen::Index num_snps = 10;
        auto [bed_prefix, expected_genotypes]
            = bed_fixture.create_bed_files(num_samples, num_snps, 0.0);

        PredictDataPipe::Config config;
        config.bed_path = bed_prefix;
        config.qcovar_path = "";
        config.dcovar_path = "";
        config.iid_only = false;

        PredictDataPipe pipe(config);

        REQUIRE(pipe.num_qcovariates() == 0);
        REQUIRE(pipe.num_dcovariates() == 0);
        REQUIRE(pipe.qcovariate_names().empty());
        REQUIRE(pipe.dcovariate_names().empty());

        auto genotypes = std::move(pipe).take_data().genotype;
        REQUIRE(genotypes.rows() == num_samples);
        REQUIRE(genotypes.cols() == num_snps);

        gelex::test::are_matrices_equal(genotypes, expected_genotypes, 1e-8);
    }
}

TEST_CASE(
    "PredictDataPipe - With quantitative covariates",
    "[predict][predict_pipe]")
{
    BedFixture bed_fixture;
    FileFixture file_fixture;

    SECTION("BED + qcovar file")
    {
        const Eigen::Index num_samples = 3;
        const Eigen::Index num_snps = 5;
        auto [bed_prefix, expected_genotypes]
            = bed_fixture.create_bed_files(num_samples, num_snps, 0.0);

        // Read FAM file to get sample IDs
        auto fam_path = bed_prefix;
        fam_path.replace_extension(".fam");
        auto [fids, iids] = read_fam(fam_path);

        // Create qcovar file
        std::string qcovar_content = "FID\tIID\tAge\tHeight\n";
        for (size_t i = 0; i < fids.size(); ++i)
        {
            qcovar_content += std::format(
                "{}\t{}\t{}\t{}\n",
                fids[i],
                iids[i],
                20 + i,
                1.6 + (static_cast<double>(i) * 0.1));
        }
        auto qcovar_path
            = file_fixture.create_text_file(qcovar_content, ".qcovar");

        PredictDataPipe::Config config;
        config.bed_path = bed_prefix;
        config.qcovar_path = qcovar_path;
        config.dcovar_path = "";
        config.iid_only = false;

        PredictDataPipe pipe(config);

        REQUIRE(pipe.num_qcovariates() == 2);
        REQUIRE(pipe.num_dcovariates() == 0);

        const auto& qnames = pipe.qcovariate_names();
        REQUIRE(qnames.size() == 2);
        REQUIRE(qnames[0] == "Age");
        REQUIRE(qnames[1] == "Height");

        auto data = std::move(pipe).take_data();
        auto& qcovariates = data.qcovariates;
        auto& genotypes = data.genotype;
        REQUIRE(qcovariates.rows() == num_samples);
        REQUIRE(qcovariates.cols() == 3);  // intercept + 2 covariates

        // Check intercept column
        for (Eigen::Index i = 0; i < qcovariates.rows(); ++i)
        {
            REQUIRE(qcovariates(i, 0) == 1.0);
        }

        // Check covariate values
        for (size_t i = 0; i < fids.size(); ++i)
        {
            REQUIRE(qcovariates(i, 1) == 20 + i);
            REQUIRE(std::abs(qcovariates(i, 2) - (1.6 + i * 0.1)) < 1e-10);
        }

        REQUIRE(genotypes.rows() == num_samples);
        REQUIRE(genotypes.cols() == num_snps);
        gelex::test::are_matrices_equal(genotypes, expected_genotypes, 1e-8);
    }
}

TEST_CASE(
    "PredictDataPipe - With categorical covariates",
    "[predict][predict_pipe]")
{
    BedFixture bed_fixture;

    SECTION("BED + covar file")
    {
        const Eigen::Index num_samples = 3;
        const Eigen::Index num_snps = 5;
        auto [bed_prefix, expected_genotypes]
            = bed_fixture.create_bed_files(num_samples, num_snps, 0.0);

        auto& file_fixture = bed_fixture.get_file_fixture();
        // Read FAM file to get sample IDs
        auto fam_path = bed_prefix;
        fam_path.replace_extension(".fam");
        auto [fids, iids] = read_fam(fam_path);

        // Create covar file
        std::string covar_content = "FID\tIID\tSex\tPopulation\n";
        for (size_t i = 0; i < fids.size(); ++i)
        {
            covar_content += std::format(
                "{}\t{}\t{}\t{}\n",
                fids[i],
                iids[i],
                (i % 2 == 0 ? "M" : "F"),
                (i % 3 == 0 ? "EUR" : "AFR"));
        }

        auto covar_path
            = file_fixture.create_text_file(covar_content, ".covar");

        PredictDataPipe::Config config;
        config.bed_path = bed_prefix;
        config.dcovar_path = covar_path;
        config.iid_only = false;

        PredictDataPipe pipe(config);

        REQUIRE(pipe.num_qcovariates() == 0);
        REQUIRE(pipe.num_dcovariates() == 2);

        const auto& cnames = pipe.dcovariate_names();
        REQUIRE(cnames.size() == 2);
        REQUIRE(cnames[0] == "Sex");
        REQUIRE(cnames[1] == "Population");

        auto data = std::move(pipe).take_data();
        auto& covariates = data.dcovariates;
        auto& genotypes = data.genotype;
        REQUIRE(covariates.size() == 2);
        REQUIRE(covariates.count("Sex") == 1);
        REQUIRE(covariates.count("Population") == 1);

        const auto& sex_values = covariates.at("Sex");
        REQUIRE(sex_values.size() == num_samples);
        for (size_t i = 0; i < sex_values.size(); ++i)
        {
            REQUIRE(sex_values[i] == (i % 2 == 0 ? "M" : "F"));
        }

        const auto& pop_values = covariates.at("Population");
        REQUIRE(pop_values.size() == num_samples);
        for (size_t i = 0; i < pop_values.size(); ++i)
        {
            REQUIRE(pop_values[i] == (i % 3 == 0 ? "EUR" : "AFR"));
        }

        REQUIRE(genotypes.rows() == num_samples);
        REQUIRE(genotypes.cols() == num_snps);
        gelex::test::are_matrices_equal(genotypes, expected_genotypes, 1e-8);
    }
}

TEST_CASE("PredictDataPipe - Sample intersection", "[predict][predict_pipe]")
{
    BedFixture bed_fixture;
    FileFixture file_fixture;

    SECTION("Partial sample overlap")
    {
        const Eigen::Index num_bed_samples = 5;
        const Eigen::Index num_snps = 5;
        auto [bed_prefix, expected_genotypes]
            = bed_fixture.create_bed_files(num_bed_samples, num_snps, 0.0);

        // Read FAM file to get sample IDs
        auto fam_path = bed_prefix;
        fam_path.replace_extension(".fam");
        auto [fids, iids] = read_fam(fam_path);

        // Create qcovar file with only first 2 samples
        std::string qcovar_content = "FID\tIID\tAge\n";
        for (size_t i = 0; i < 2; ++i)
        {
            qcovar_content += fids[i] + "\t" + iids[i] + "\t";
            qcovar_content += std::to_string(20 + i) + "\n";
        }

        auto qcovar_path
            = file_fixture.create_text_file(qcovar_content, ".qcovar");

        PredictDataPipe::Config config;
        config.bed_path = bed_prefix;
        config.qcovar_path = qcovar_path;
        config.dcovar_path = "";
        config.iid_only = false;

        PredictDataPipe pipe(config);

        REQUIRE(pipe.num_qcovariates() == 1);
        auto data = std::move(pipe).take_data();
        auto& qcovariates = data.qcovariates;
        auto& genotypes = data.genotype;
        REQUIRE(qcovariates.rows() == 2);  // Only 2 common samples
        REQUIRE(qcovariates.cols() == 2);  // intercept + 1 covariate
        REQUIRE(qcovariates(0, 0) == 1.0);  // intercept
        REQUIRE(qcovariates(1, 0) == 1.0);  // intercept
        REQUIRE(qcovariates(0, 1) == 20);   // Age for first sample
        REQUIRE(qcovariates(1, 1) == 21);   // Age for second sample

        REQUIRE(genotypes.rows() == 2);
        REQUIRE(genotypes.cols() == num_snps);

        // Check that genotypes are for first 2 samples
        Eigen::MatrixXd expected_subset = expected_genotypes.topRows(2);
        gelex::test::are_matrices_equal(genotypes, expected_subset, 1e-8);
    }
}

TEST_CASE("PredictDataPipe - iid_only mode", "[predict][predict_pipe]")
{
    BedFixture bed_fixture;
    FileFixture file_fixture;

    SECTION("iid_only = true")
    {
        const Eigen::Index num_samples = 3;
        const Eigen::Index num_snps = 5;
        auto [bed_prefix, expected_genotypes]
            = bed_fixture.create_bed_files(num_samples, num_snps, 0.0);

        // Read FAM file to get IIDs
        auto fam_path = bed_prefix;
        fam_path.replace_extension(".fam");
        auto [_, iids] = read_fam(fam_path);

        // Create qcovar file using only IIDs
        std::string qcovar_content = "FID\tIID\tAge\n";
        for (size_t i = 0; i < iids.size(); ++i)
        {
            qcovar_content += "1\t" + iids[i] + "\t";
            qcovar_content += std::to_string(20 + i) + "\n";
        }

        auto qcovar_path
            = file_fixture.create_text_file(qcovar_content, ".qcovar");

        PredictDataPipe::Config config;
        config.bed_path = bed_prefix;
        config.qcovar_path = qcovar_path;
        config.dcovar_path = "";
        config.iid_only = true;

        PredictDataPipe pipe(config);

        REQUIRE(pipe.num_qcovariates() == 1);

        auto data = std::move(pipe).take_data();
        auto& qcovariates = data.qcovariates;
        auto& genotypes = data.genotype;
        REQUIRE(qcovariates.rows() == num_samples);
        REQUIRE(qcovariates.cols() == 2);  // intercept + Age

        // Check intercept column
        for (Eigen::Index i = 0; i < qcovariates.rows(); ++i)
        {
            REQUIRE(qcovariates(i, 0) == 1.0);
        }

        // Check Age values
        for (size_t i = 0; i < iids.size(); ++i)
        {
            REQUIRE(qcovariates(i, 1) == 20 + i);
        }

        REQUIRE(genotypes.rows() == num_samples);
        REQUIRE(genotypes.cols() == num_snps);
        gelex::test::are_matrices_equal(genotypes, expected_genotypes, 1e-8);
    }
}

TEST_CASE(
    "PredictDataPipe - Data movement methods",
    "[predict][predict_pipe]")
{
    BedFixture bed_fixture;
    FileFixture file_fixture;

    SECTION("take methods move data out")
    {
        const Eigen::Index num_samples = 3;
        const Eigen::Index num_snps = 5;
        auto [bed_prefix, expected_genotypes]
            = bed_fixture.create_bed_files(num_samples, num_snps, 0.0);

        // Read FAM file
        auto fam_path = bed_prefix;
        fam_path.replace_extension(".fam");
        auto [fids, iids] = read_fam(fam_path);

        // Create qcovar file
        std::string qcovar_content = "FID\tIID\tAge\n";
        for (size_t i = 0; i < fids.size(); ++i)
        {
            qcovar_content += fids[i] + "\t" + iids[i] + "\t";
            qcovar_content += std::to_string(20 + i) + "\n";
        }

        auto qcovar_path
            = file_fixture.create_text_file(qcovar_content, ".qcovar");

        PredictDataPipe::Config config;
        config.bed_path = bed_prefix;
        config.qcovar_path = qcovar_path;
        config.dcovar_path = "";
        config.iid_only = false;

        PredictDataPipe pipe(config);

        // Take all three components
        auto data = std::move(pipe).take_data();
        auto& qcovariates = data.qcovariates;
        auto& genotypes = data.genotype;

        REQUIRE(qcovariates.rows() == num_samples);
        REQUIRE(qcovariates.cols() == 2);
        REQUIRE(genotypes.rows() == num_samples);
        REQUIRE(genotypes.cols() == num_snps);

        gelex::test::are_matrices_equal(genotypes, expected_genotypes, 1e-8);
    }
}

TEST_CASE("PredictDataPipe - Edge cases", "[predict][predict_pipe]")
{
    BedFixture bed_fixture;
    FileFixture file_fixture;

    SECTION("Empty qcovar file")
    {
        const Eigen::Index num_samples = 3;
        const Eigen::Index num_snps = 5;
        auto [bed_prefix, expected_genotypes]
            = bed_fixture.create_bed_files(num_samples, num_snps, 0.0);

        // Read FAM file to get sample IDs
        auto fam_path = bed_prefix;
        fam_path.replace_extension(".fam");
        auto [fids, iids] = read_fam(fam_path);

        // Create qcovar file with one dummy covariate column
        std::string qcovar_content = "FID\tIID\tDummy\n";
        for (size_t i = 0; i < fids.size(); ++i)
        {
            qcovar_content += fids[i] + "\t" + iids[i] + "\t0.0\n";
        }

        auto qcovar_path
            = file_fixture.create_text_file(qcovar_content, ".qcovar");

        PredictDataPipe::Config config;
        config.bed_path = bed_prefix;
        config.qcovar_path = qcovar_path;
        config.dcovar_path = "";
        config.iid_only = false;

        PredictDataPipe pipe(config);

        REQUIRE(pipe.num_qcovariates() == 1);

        auto data = std::move(pipe).take_data();
        auto& qcovariates = data.qcovariates;
        auto& genotypes = data.genotype;
        REQUIRE(qcovariates.rows() == num_samples);
        REQUIRE(qcovariates.cols() == 2);  // intercept + dummy

        // Check intercept column
        for (Eigen::Index i = 0; i < qcovariates.rows(); ++i)
        {
            REQUIRE(qcovariates(i, 0) == 1.0);
        }
        // Check dummy column (all zeros)
        for (Eigen::Index i = 0; i < qcovariates.rows(); ++i)
        {
            REQUIRE(qcovariates(i, 1) == 0.0);
        }

        REQUIRE(genotypes.rows() == num_samples);
        REQUIRE(genotypes.cols() == num_snps);
        gelex::test::are_matrices_equal(genotypes, expected_genotypes, 1e-8);
    }

    SECTION("No common samples")
    {
        const Eigen::Index num_samples = 3;
        const Eigen::Index num_snps = 5;
        auto [bed_prefix, expected_genotypes]
            = bed_fixture.create_bed_files(num_samples, num_snps, 0.0);

        // Create qcovar file with non-existent samples
        auto qcovar_path = file_fixture.create_text_file(
            "FID\tIID\tAge\n"
            "999\t999\t20\n"
            "888\t888\t21\n",
            ".qcovar");

        PredictDataPipe::Config config;
        config.bed_path = bed_prefix;
        config.qcovar_path = qcovar_path;
        config.dcovar_path = "";
        config.iid_only = false;

        PredictDataPipe pipe(config);

        REQUIRE(pipe.num_qcovariates() == 1);
        auto data = std::move(pipe).take_data();
        auto& qcovariates = data.qcovariates;
        auto& genotypes = data.genotype;
        REQUIRE(qcovariates.rows() == 0);  // No common samples
        REQUIRE(qcovariates.cols() == 2);

        REQUIRE(genotypes.rows() == 0);
        REQUIRE(genotypes.cols() == num_snps);
    }
}
