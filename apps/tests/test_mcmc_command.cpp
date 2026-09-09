// Copyright 2026 RuLei Chen
// SPDX-License-Identifier: Apache-2.0

#include <Eigen/Core>
#include <catch2/catch_test_macros.hpp>
#include <cstdint>
#include <filesystem>
#include <string>

#include "gelex/bayes/builtin_method.h"
#include "gelex/bayes/draws.h"
#include "gelex/data/genotype_method.h"
#include "gelex/data/snp_lut_io.h"
#include "gelex/exception.h"
#include "gelex/genetic_mode.h"
#include "gelex/io/binary_format.h"
#include "gelex/io/csc_reader.h"
#include "gelex/io/dense_reader.h"

#include "bed_fixture.h"
#include "cli/logging.h"
#include "cli/mcmc/command.h"
#include "cli/mcmc/config.h"
#include "cli/runtime.h"

namespace cli
{

auto setup_parallelization(int /*num_threads*/) -> void {}

}  // namespace cli

TEST_CASE(
    "MCMC command runs typed fitting with discrete and quantitative random "
    "effects",
    "[cli][mcmc][command]")
{
    gelex::test::BedFixture fixture;
    const auto [bfile, genotypes] = fixture.create_deterministic_bed_files(
        Eigen::MatrixXd{
            {0.0, 0.0, 0.0},
            {1.0, 0.0, 1.0},
            {2.0, 1.0, 0.0},
            {0.0, 2.0, 1.0},
            {2.0, 2.0, 2.0}},
        {"I1", "I2", "I3", "I4", "I5"});
    static_cast<void>(genotypes);
    auto& files = fixture.get_file_fixture();
    const auto phenotype = files.create_named_text_file(
        "phenotype.tsv",
        "FID\tIID\tTrait\n"
        "fam3\tI3\t2.0\n"
        "fam1\tI1\t1.0\n"
        "fam5\tI5\t4.0\n"
        "fam2\tI2\t3.0\n"
        "fam4\tI4\t5.0\n");
    const auto drand = files.create_named_text_file(
        "random_factors.tsv",
        "FID\tIID\tGroup\n"
        "fam2\tI2\tB\n"
        "fam4\tI4\tB\n"
        "fam1\tI1\tA\n"
        "fam5\tI5\tA\n"
        "fam3\tI3\tA\n");
    const auto qrand = files.create_named_text_file(
        "random_slopes.tsv",
        "FID\tIID\tSlope\n"
        "fam5\tI5\t2.0\n"
        "fam2\tI2\t-1.0\n"
        "fam4\tI4\t1.0\n"
        "fam1\tI1\t-2.0\n"
        "fam3\tI3\t0.0\n");
    const auto output = files.get_test_dir() / "typed_mcmc";

    cli::McmcConfig config;
    config.base_data.pheno_path = phenotype.string();
    config.random.drand_path = drand.string();
    config.random.qrand_paths = {qrand.string()};
    config.random_pve = 0.1;
    config.bfile = bfile.string();
    config.out = output.string();
    config.iters = 4;
    config.burn_in = 2;
    config.thin = 1;
    config.threads = 1;
    cli::logging::initialize(config.out);

    REQUIRE(mcmc_execute(config) == 0);

    const gelex::DenseReader draws(config.out + ".draws");
    REQUIRE(draws.contains("random/Group/coefficients"));
    REQUIRE(draws.contains("random/random_slopes/coefficients"));
    REQUIRE(draws.to_map<float>("random/Group/coefficients").rows() == 2);
    REQUIRE(draws.to_map<float>("random/Group/coefficients").cols() == 2);
    REQUIRE(
        draws.to_map<float>("random/random_slopes/coefficients").rows() == 1);
    REQUIRE(
        draws.to_map<float>("random/random_slopes/coefficients").cols() == 2);

    REQUIRE(draws.to_map<double>("fixed/coefficients").cols() == 2);
    REQUIRE(draws.to_map<double>("random/Group/variance").allFinite());
    REQUIRE(draws.to_map<double>("random/random_slopes/variance").allFinite());
    REQUIRE(draws.to_map<float>("genetic/A/coefficients").rows() == 3);
    REQUIRE(draws.to_map<float>("genetic/A/coefficients").cols() == 2);

    const auto luts = gelex::load_snp_luts(config.out + ".snplut");
    REQUIRE(luts.size() == 1);
    REQUIRE(luts.contains(gelex::GeneticMode::A));
    REQUIRE(luts.at(gelex::GeneticMode::A).cols() == 3);
}

TEST_CASE(
    "MCMC command fits BayesCD with marker annotations",
    "[cli][mcmc][command]")
{
    gelex::test::BedFixture fixture;
    const auto [bfile, genotypes] = fixture.create_deterministic_bed_files(
        Eigen::MatrixXd{
            {0.0, 0.0, 1.0},
            {1.0, 0.0, 2.0},
            {2.0, 1.0, 0.0},
            {0.0, 2.0, 1.0},
            {1.0, 2.0, 0.0},
            {2.0, 1.0, 2.0}},
        {"I1", "I2", "I3", "I4", "I5", "I6"},
        {"snp1", "snp2", "snp3"},
        {"1", "1", "1"},
        {{'A', 'G'}, {'A', 'G'}, {'A', 'G'}});
    static_cast<void>(genotypes);
    auto& files = fixture.get_file_fixture();
    const auto phenotype = files.create_named_text_file(
        "phenotype.tsv",
        "FID\tIID\tTrait\n"
        "fam1\tI1\t1.0\n"
        "fam2\tI2\t-0.5\n"
        "fam3\tI3\t0.25\n"
        "fam4\tI4\t2.0\n"
        "fam5\tI5\t-1.0\n"
        "fam6\tI6\t0.75\n");
    const auto output = files.get_test_dir() / "bayescd";

    cli::McmcConfig config;
    config.base_data.pheno_path = phenotype.string();
    config.bfile = bfile.string();
    config.out = output.string();
    config.mode = gelex::GeneticMode::A | gelex::GeneticMode::D;
    config.method = gelex::BayesMethod::CD;
    config.geno_method = gelex::GenotypeMethod::NOIACenter;
    config.iters = 4;
    config.burn_in = 2;
    config.thin = 1;
    config.threads = 1;
    cli::logging::initialize(config.out);

    SECTION("single annotation column")
    {
        config.manno = files
                           .create_named_text_file(
                               "markers.anno",
                               "CHR\tSNP\tBP\tA1\tA2\tAnnotation\n"
                               "1\tsnp1\t1\tA\tG\t-0.5\n"
                               "1\tsnp2\t2\tA\tG\t0.0\n"
                               "1\tsnp3\t3\tA\tG\t0.5\n")
                           .string();

        REQUIRE(mcmc_execute(config) == 0);

        const auto draws_path = config.out + ".draws";
        const gelex::DenseReader draws(draws_path);
        const gelex::CscReader sparse(gelex::sparse_draws_path(draws_path));
        REQUIRE(sparse.contains("genetic/A/coefficients"));
        REQUIRE(sparse.contains("genetic/D/coefficients"));
        REQUIRE(draws.contains("genetic/D/annotation_coefficients"));
        REQUIRE(sparse.contains("genetic/joint/assignment"));

        REQUIRE(
            draws.to_map<float>("genetic/D/annotation_coefficients").rows()
            == 2);
        REQUIRE(
            draws.to_map<float>("genetic/D/annotation_coefficients").cols()
            == 2);
        REQUIRE(
            sparse.info("genetic/joint/assignment").shape
            == (gelex::BinaryShape{3, 2}));

        const auto luts = gelex::load_snp_luts(config.out + ".snplut");
        REQUIRE(luts.size() == 2);
        REQUIRE(luts.contains(gelex::GeneticMode::A));
        REQUIRE(luts.contains(gelex::GeneticMode::D));
    }

    SECTION("multiple annotation columns are rejected")
    {
        config.manno = files
                           .create_named_text_file(
                               "markers.anno",
                               "CHR\tSNP\tBP\tA1\tA2\tFirst\tSecond\n"
                               "1\tsnp1\t1\tA\tG\t-0.5\t1.0\n"
                               "1\tsnp2\t2\tA\tG\t0.0\t1.0\n"
                               "1\tsnp3\t3\tA\tG\t0.5\t1.0\n")
                           .string();
        REQUIRE_THROWS_AS(mcmc_execute(config), gelex::GelexException);
    }
}
