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

#include <algorithm>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <random>
#include <sstream>
#include <string>
#include <tuple>
#include <unordered_map>

#include <Eigen/Core>
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include "bed_fixture.h"
#include "gelex/data/genotype_processor.h"
#include "gelex/data/simulate.h"
#include "gelex/exception.h"
#include "utils/math_utils.h"

namespace fs = std::filesystem;

using namespace gelex;  // NOLINT
using Catch::Matchers::WithinAbs;
using gelex::test::BedFixture;

namespace
{

constexpr double VARIANCE_TOLERANCE = 0.1;

struct CausalEffect
{
    double additive = 0.0;
    double dominance = 0.0;
    int add_class = 0;
    int dom_class = 0;
};

auto count_lines(const fs::path& path) -> int
{
    std::ifstream file(path);
    return static_cast<int>(std::count(
        std::istreambuf_iterator<char>(file),
        std::istreambuf_iterator<char>(),
        '\n'));
}

auto read_first_line(const fs::path& path) -> std::string
{
    std::ifstream file(path);
    std::string line;
    std::getline(file, line);
    return line;
}

auto read_file_content(const fs::path& path) -> std::string
{
    std::ifstream file(path);
    std::ostringstream ss;
    ss << file.rdbuf();
    return ss.str();
}

auto make_config(
    const fs::path& bed_path,
    double h2 = 0.5,
    double d2 = 0.0,
    int seed = 42,
    std::vector<gelex::EffectSizeClass> add_classes = {{1.0, 1.0}},
    std::vector<gelex::EffectSizeClass> dom_classes = {{1.0, 1.0}})
    -> PhenotypeSimulator::Config
{
    return {
        .bed_path = bed_path,
        .add_heritability = h2,
        .dom_heritability = d2,
        .add_effect_classes = std::move(add_classes),
        .dom_effect_classes = std::move(dom_classes),
        .seed = seed,
        .output_path = {},
    };
}

auto parse_causal_effects(const fs::path& causal_path)
    -> std::unordered_map<std::string, CausalEffect>
{
    std::ifstream file(causal_path);
    std::string line;
    std::getline(file, line);  // skip header

    std::unordered_map<std::string, CausalEffect> effects;
    while (std::getline(file, line))
    {
        std::istringstream iss(line);
        std::string snp;
        CausalEffect effect;
        iss >> snp >> effect.additive >> effect.dominance >> effect.add_class
            >> effect.dom_class;
        effects[snp] = effect;
    }
    return effects;
}

auto parse_snp_indices(const fs::path& bim_path)
    -> std::unordered_map<std::string, Eigen::Index>
{
    std::ifstream file(bim_path);
    std::string line;
    std::unordered_map<std::string, Eigen::Index> snp_to_col;
    Eigen::Index col = 0;

    while (std::getline(file, line))
    {
        std::istringstream iss(line);
        std::string chr;
        std::string snp_id;
        iss >> chr >> snp_id;
        snp_to_col[snp_id] = col++;
    }
    return snp_to_col;
}

auto parse_phenotypes(const fs::path& phen_path, Eigen::Index n_samples)
    -> Eigen::VectorXd
{
    std::ifstream file(phen_path);
    std::string line;
    std::getline(file, line);  // skip header

    Eigen::VectorXd phenotypes(n_samples);
    Eigen::Index idx = 0;
    while (std::getline(file, line))
    {
        std::istringstream iss(line);
        std::string fid;
        std::string iid;
        double val{};
        iss >> fid >> iid >> val;
        phenotypes(idx++) = val;
    }
    return phenotypes;
}

auto generate_random_genotypes(
    Eigen::Index n_samples,
    Eigen::Index n_snps,
    unsigned int seed) -> Eigen::MatrixXd
{
    std::mt19937_64 rng(seed);
    Eigen::MatrixXd genotypes(n_samples, n_snps);

    for (Eigen::Index j = 0; j < n_snps; ++j)
    {
        double maf
            = 0.1
              + (0.3 * static_cast<double>(j) / static_cast<double>(n_snps));
        std::binomial_distribution<int> binom(2, maf);
        for (Eigen::Index i = 0; i < n_samples; ++i)
        {
            genotypes(i, j) = static_cast<double>(binom(rng));
        }
    }
    return genotypes;
}

auto extract_causal_columns(
    const Eigen::MatrixXd& genotypes,
    const std::unordered_map<std::string, CausalEffect>& effects,
    const std::unordered_map<std::string, Eigen::Index>& snp_to_col)
    -> std::tuple<Eigen::MatrixXd, Eigen::VectorXd, Eigen::VectorXd>
{
    auto n_causal = static_cast<Eigen::Index>(effects.size());
    Eigen::MatrixXd causal_geno(genotypes.rows(), n_causal);
    Eigen::VectorXd add_betas(n_causal);
    Eigen::VectorXd dom_betas(n_causal);

    Eigen::Index i = 0;
    for (const auto& [snp, effect] : effects)
    {
        auto it = snp_to_col.find(snp);
        REQUIRE(it != snp_to_col.end());
        causal_geno.col(i) = genotypes.col(it->second);
        add_betas(i) = effect.additive;
        dom_betas(i) = effect.dominance;
        ++i;
    }
    return {causal_geno, add_betas, dom_betas};
}

}  // namespace

TEST_CASE("PhenotypeSimulator - parameter validation", "[simulate]")
{
    BedFixture fixture;
    auto [bed_path, _] = fixture.create_bed_files(10, 20, 0.0, 0.05, 0.5, 42);

    SECTION("Valid config does not throw")
    {
        REQUIRE_NOTHROW(PhenotypeSimulator(make_config(bed_path)));
    }

    SECTION("h2 must be in (0, 1)")
    {
        REQUIRE_THROWS_AS(
            PhenotypeSimulator(make_config(bed_path, 0.0)),
            ArgumentValidationException);
        REQUIRE_THROWS_AS(
            PhenotypeSimulator(make_config(bed_path, -0.1)),
            ArgumentValidationException);
        REQUIRE_THROWS_AS(
            PhenotypeSimulator(make_config(bed_path, 1.0)),
            ArgumentValidationException);
        REQUIRE_THROWS_AS(
            PhenotypeSimulator(make_config(bed_path, 1.5)),
            ArgumentValidationException);
    }

    SECTION("d2 must be in [0, 1)")
    {
        REQUIRE_THROWS_AS(
            PhenotypeSimulator(make_config(bed_path, 0.5, -0.1)),
            ArgumentValidationException);
        REQUIRE_THROWS_AS(
            PhenotypeSimulator(make_config(bed_path, 0.5, 1.0)),
            ArgumentValidationException);
    }

    SECTION("h2 + d2 must be < 1")
    {
        REQUIRE_THROWS_AS(
            PhenotypeSimulator(make_config(bed_path, 0.6, 0.5)),
            ArgumentValidationException);
    }
}

TEST_CASE("PhenotypeSimulator - basic simulation", "[simulate]")
{
    BedFixture fixture;
    constexpr Eigen::Index N_SAMPLES = 50;
    constexpr Eigen::Index N_SNPS = 100;
    auto [bed_path, _]
        = fixture.create_bed_files(N_SAMPLES, N_SNPS, 0.0, 0.05, 0.5, 42);

    SECTION("Default output generates .phen and .causal files")
    {
        PhenotypeSimulator(make_config(bed_path)).simulate();

        REQUIRE(fs::exists(fs::path(bed_path).replace_extension(".phen")));
        REQUIRE(fs::exists(fs::path(bed_path).replace_extension(".causal")));
    }

    SECTION("Custom output path")
    {
        auto output_path
            = fixture.get_file_fixture().get_test_dir() / "custom_output.phen";
        auto config = make_config(bed_path);
        config.output_path = output_path;

        PhenotypeSimulator(config).simulate();

        REQUIRE(fs::exists(output_path));
        REQUIRE(fs::exists(fs::path(output_path).replace_extension(".causal")));
    }

    SECTION("Phenotype file format")
    {
        PhenotypeSimulator(make_config(bed_path)).simulate();

        auto phen_path = fs::path(bed_path).replace_extension(".phen");
        REQUIRE(read_first_line(phen_path) == "FID\tIID\tphenotype");
        REQUIRE(count_lines(phen_path) == N_SAMPLES + 1);  // header + samples
    }

    SECTION("Causal file format")
    {
        PhenotypeSimulator(make_config(bed_path)).simulate();

        auto causal_path = fs::path(bed_path).replace_extension(".causal");
        REQUIRE(
            read_first_line(causal_path)
            == "SNP\tadditive_effect\tdominance_effect\t"
               "add_class\tdom_class");

        REQUIRE(count_lines(causal_path) == N_SNPS + 1);
    }
}

TEST_CASE("PhenotypeSimulator - reproducibility", "[simulate]")
{
    BedFixture fixture;
    auto [bed_path, _] = fixture.create_bed_files(50, 100, 0.0, 0.05, 0.5, 42);

    auto output1
        = fixture.get_file_fixture().get_test_dir() / "repro_run1.phen";
    auto output2
        = fixture.get_file_fixture().get_test_dir() / "repro_run2.phen";

    auto config = make_config(bed_path, 0.5, 0.0, 123);
    config.output_path = output1;
    PhenotypeSimulator(config).simulate();

    config.output_path = output2;
    PhenotypeSimulator(config).simulate();

    REQUIRE(read_file_content(output1) == read_file_content(output2));
    REQUIRE(
        read_file_content(fs::path(output1).replace_extension(".causal"))
        == read_file_content(fs::path(output2).replace_extension(".causal")));
}

TEST_CASE("PhenotypeSimulator - dominance effects", "[simulate]")
{
    BedFixture fixture;
    auto [bed_path, _] = fixture.create_bed_files(50, 100, 0.0, 0.05, 0.5, 42);

    PhenotypeSimulator(make_config(bed_path, 0.5, 0.2)).simulate();

    auto causal_path = fs::path(bed_path).replace_extension(".causal");
    REQUIRE(fs::exists(causal_path));

    auto effects = parse_causal_effects(causal_path);
    bool has_nonzero_dominance = std::any_of(
        effects.begin(),
        effects.end(),
        [](const auto& kv) { return std::abs(kv.second.dominance) > 1e-10; });

    REQUIRE(has_nonzero_dominance);
}

TEST_CASE("PhenotypeSimulator - additive variance", "[simulate]")
{
    BedFixture fixture;
    constexpr Eigen::Index N_SAMPLES = 500;
    constexpr Eigen::Index N_SNPS = 200;
    constexpr double H2 = 0.5;

    auto genotypes = generate_random_genotypes(N_SAMPLES, N_SNPS, 99);
    auto [bed_path, stored_geno]
        = fixture.create_deterministic_bed_files(genotypes);

    PhenotypeSimulator(make_config(bed_path, H2)).simulate();

    auto effects
        = parse_causal_effects(fs::path(bed_path).replace_extension(".causal"));
    auto snp_to_col
        = parse_snp_indices(fs::path(bed_path).replace_extension(".bim"));
    auto [causal_geno, add_betas, _]
        = extract_causal_columns(stored_geno, effects, snp_to_col);

    process_matrix<grm::OrthStandardized::Additive>(causal_geno);
    Eigen::VectorXd g_a = causal_geno * add_betas;

    auto phenotypes = parse_phenotypes(
        fs::path(bed_path).replace_extension(".phen"), N_SAMPLES);

    double observed_h2 = detail::var(g_a)(0) / detail::var(phenotypes)(0);
    REQUIRE_THAT(observed_h2, WithinAbs(H2, VARIANCE_TOLERANCE));
}

TEST_CASE("PhenotypeSimulator - additive and dominance variance", "[simulate]")
{
    BedFixture fixture;
    constexpr Eigen::Index N_SAMPLES = 500;
    constexpr Eigen::Index N_SNPS = 200;
    constexpr double H2 = 0.4;
    constexpr double D2 = 0.2;

    auto genotypes = generate_random_genotypes(N_SAMPLES, N_SNPS, 99);
    auto [bed_path, stored_geno]
        = fixture.create_deterministic_bed_files(genotypes);

    PhenotypeSimulator(make_config(bed_path, H2, D2)).simulate();

    auto effects
        = parse_causal_effects(fs::path(bed_path).replace_extension(".causal"));
    auto snp_to_col
        = parse_snp_indices(fs::path(bed_path).replace_extension(".bim"));
    auto [causal_geno, add_betas, dom_betas]
        = extract_causal_columns(stored_geno, effects, snp_to_col);

    // Additive genetic values
    Eigen::MatrixXd x_add = causal_geno;
    process_matrix<grm::OrthStandardized::Additive>(x_add);
    Eigen::VectorXd g_a = x_add * add_betas;

    // Dominance genetic values
    Eigen::MatrixXd x_dom = causal_geno;
    process_matrix<grm::OrthStandardized::Dominant>(x_dom);
    Eigen::VectorXd g_d = x_dom * dom_betas;

    // Scale dominance: scaled_d = d * sqrt(target / raw), target = Va * d2 / h2
    double var_ga = detail::var(g_a)(0);
    double var_gd_raw = detail::var(g_d)(0);
    double scale = std::sqrt(var_ga * D2 / H2 / var_gd_raw);
    Eigen::VectorXd g_d_scaled = g_d * scale;

    auto phenotypes = parse_phenotypes(
        fs::path(bed_path).replace_extension(".phen"), N_SAMPLES);

    double var_phen = detail::var(phenotypes)(0);
    REQUIRE_THAT(var_ga / var_phen, WithinAbs(H2, VARIANCE_TOLERANCE));
    REQUIRE_THAT(
        detail::var(g_d_scaled)(0) / var_phen,
        WithinAbs(D2, VARIANCE_TOLERANCE));
}

TEST_CASE("PhenotypeSimulator - mixture normal effect classes", "[simulate]")
{
    BedFixture fixture;
    constexpr Eigen::Index N_SAMPLES = 200;
    constexpr Eigen::Index N_SNPS = 200;

    auto [bed_path, _]
        = fixture.create_bed_files(N_SAMPLES, N_SNPS, 0.0, 0.05, 0.5, 42);

    // 3-class mixture: small/medium/large effect sizes
    std::vector<gelex::EffectSizeClass> add_classes
        = {{0.5, 0.0001}, {0.3, 0.01}, {0.2, 1.0}};

    auto config = make_config(bed_path, 0.5, 0.0, 42, add_classes);
    PhenotypeSimulator(config).simulate();

    auto causal_path = fs::path(bed_path).replace_extension(".causal");
    auto effects = parse_causal_effects(causal_path);

    REQUIRE(static_cast<int>(effects.size()) == N_SNPS);

    // Count SNPs per class
    std::array<int, 3> class_counts = {0, 0, 0};
    for (const auto& [snp, effect] : effects)
    {
        REQUIRE(effect.add_class >= 0);
        REQUIRE(effect.add_class < 3);
        class_counts[effect.add_class]++;
    }

    // Each class should have at least one SNP
    for (int cls = 0; cls < 3; ++cls)
    {
        REQUIRE(class_counts[cls] > 0);
    }

    // Class 0 (50%) should have more SNPs than class 2 (20%)
    REQUIRE(class_counts[0] > class_counts[2]);
}

// Effect class validation is now handled by EffectSampler (see test_effect_sampler.cpp)
