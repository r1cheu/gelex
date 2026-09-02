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

#include <Eigen/Core>
#include <catch2/catch_test_macros.hpp>
#include <filesystem>
#include <fstream>
#include <iterator>
#include <string>

#include "gelex/data/reader.h"
#include "gelex/data/snp_lut_io.h"
#include "gelex/genetic_mode.h"
#include "gelex/io/binary_reader.h"

#include "bed_fixture.h"
#include "cli/logging.h"
#include "cli/mcmc/command.h"
#include "cli/mcmc/config.h"
#include "cli/runtime.h"

namespace cli
{

auto setup_parallelization(int /*num_threads*/) -> void {}

}  // namespace cli

namespace
{

auto read_text(const std::filesystem::path& path) -> std::string
{
    std::ifstream input(path);
    return {
        std::istreambuf_iterator<char>{input},
        std::istreambuf_iterator<char>{}};
}

}  // namespace

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

    const gelex::BinaryReader draws(config.out + ".draws");
    REQUIRE(draws.contains("random/Group/coefficients"));
    REQUIRE(draws.contains("random/random_slopes/coefficients"));
    REQUIRE(draws.to_map<float>("random/Group/coefficients").rows() == 2);
    REQUIRE(draws.to_map<float>("random/Group/coefficients").cols() == 2);
    REQUIRE(
        draws.to_map<float>("random/random_slopes/coefficients").rows() == 1);
    REQUIRE(
        draws.to_map<float>("random/random_slopes/coefficients").cols() == 2);

    const auto params = read_text(config.out + ".params");
    REQUIRE(params.contains("fixed/coefficients\t0\tIntercept"));
    REQUIRE(params.contains("random/Group/coefficients\t0\tGroup"));
    REQUIRE(params.contains("random/random_slopes/coefficients\t0\tSlope"));

    const auto summary = read_text(config.out + ".summary");
    REQUIRE(summary.contains("random/Group/variance\t0"));
    REQUIRE(summary.contains("random/random_slopes/variance\t0"));

    const auto snp_effects = gelex::read_snp_effects(config.out + ".snpeff");
    REQUIRE(snp_effects.rows() == 3);
    REQUIRE(snp_effects.contains("BETA_A"));
    REQUIRE(snp_effects.contains("SE_A"));
    REQUIRE(snp_effects.contains("PVE_A"));
    REQUIRE_FALSE(snp_effects.contains("PIP_A"));

    const auto luts = gelex::load_snp_luts(config.out + ".snplut");
    REQUIRE(luts.size() == 1);
    REQUIRE(luts.contains(gelex::GeneticMode::A));
    REQUIRE(luts.at(gelex::GeneticMode::A).cols() == 3);
}
