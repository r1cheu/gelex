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
#include <optional>
#include <random>
#include <sstream>
#include <string>
#include <tuple>
#include <unordered_map>
#include <utility>
#include <vector>

#include <Eigen/Core>
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include "bed_fixture.h"
#include "gelex/data/genotype/genotype_processor.h"
#include "gelex/infra/stats/descriptive.h"
#include "gelex/pipeline/simulation_engine.h"
#include "gelex/types/genotype_process_method.h"

namespace fs = std::filesystem;

using namespace gelex;  // NOLINT
using Catch::Matchers::WithinAbs;
using gelex::test::BedFixture;

namespace
{

constexpr double VARIANCE_TOLERANCE = 0.1;

struct CausalRow
{
    double additive = 0.0;
    double dominance = 0.0;
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
    Eigen::Index n_snps,
    double h2 = 0.5,
    double d2 = 0.0,
    int seed = 42,
    std::vector<EffectSize> add_classes = {},
    std::vector<EffectSize> dom_classes = {}) -> SimulationEngine::Config
{
    if (add_classes.empty())
    {
        add_classes = {{n_snps, 1.0}};
    }
    if (dom_classes.empty())
    {
        dom_classes = {{n_snps, 1.0}};
    }

    std::optional<SimulationEngine::SimulateScheme> dominance;
    if (d2 > 0.0)
    {
        dominance = SimulationEngine::SimulateScheme{
            .heritability = d2,
            .effect_sizes = std::move(dom_classes),
        };
    }
    return {
        .seed = seed,
        .bfile_prefix = bed_path,
        .output_prefix = bed_path,
        .geno_method = GenotypeProcessMethod::OrthStandardize(),
        .additive
        = SimulationEngine::SimulateScheme{
            .heritability = h2,
            .effect_sizes = std::move(add_classes),
        },
        .dominance = std::move(dominance),
    };
}

auto parse_causal_effects(const fs::path& causal_path, bool has_dominance)
    -> std::unordered_map<std::string, CausalRow>
{
    std::ifstream file(causal_path);
    std::string line;
    std::getline(file, line);

    std::unordered_map<std::string, CausalRow> effects;
    while (std::getline(file, line))
    {
        std::istringstream iss(line);
        std::string snp;
        CausalRow effect;
        iss >> snp >> effect.additive;
        if (has_dominance)
        {
            iss >> effect.dominance;
        }
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
    std::getline(file, line);

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
    const std::unordered_map<std::string, CausalRow>& effects,
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

TEST_CASE("SimulationEngine - basic output files", "[simulate]")
{
    BedFixture fixture;
    constexpr Eigen::Index N_SAMPLES = 50;
    constexpr Eigen::Index N_SNPS = 100;
    auto [bed_path, _]
        = fixture.create_bed_files(N_SAMPLES, N_SNPS, 0.0, 0.05, 0.5, 42);

    SECTION("Generates .phen and .causal files")
    {
        SimulationEngine(make_config(bed_path, N_SNPS)).run();

        REQUIRE(fs::exists(fs::path(bed_path).replace_extension(".phen")));
        REQUIRE(fs::exists(fs::path(bed_path).replace_extension(".causal")));
    }

    SECTION("Phenotype file format")
    {
        SimulationEngine(make_config(bed_path, N_SNPS)).run();

        auto phen_path = fs::path(bed_path).replace_extension(".phen");
        REQUIRE(read_first_line(phen_path) == "FID\tIID\tPhenotype");
        REQUIRE(count_lines(phen_path) == N_SAMPLES + 1);
    }

    SECTION("Causal file format (additive only)")
    {
        SimulationEngine(make_config(bed_path, N_SNPS)).run();

        auto causal_path = fs::path(bed_path).replace_extension(".causal");
        REQUIRE(read_first_line(causal_path) == "id\tadditive");
        REQUIRE(count_lines(causal_path) == N_SNPS + 1);
    }

    SECTION("Causal file format (additive + dominance)")
    {
        SimulationEngine(make_config(bed_path, N_SNPS, 0.4, 0.2)).run();

        auto causal_path = fs::path(bed_path).replace_extension(".causal");
        REQUIRE(read_first_line(causal_path) == "id\tadditive\tdominance");
        REQUIRE(count_lines(causal_path) == N_SNPS + 1);
    }
}

TEST_CASE("SimulationEngine - reproducibility", "[simulate]")
{
    BedFixture fixture1;
    BedFixture fixture2;
    constexpr Eigen::Index N_SNPS = 100;
    auto [bed_path1, _1]
        = fixture1.create_bed_files(50, N_SNPS, 0.0, 0.05, 0.5, 42);
    auto [bed_path2, _2]
        = fixture2.create_bed_files(50, N_SNPS, 0.0, 0.05, 0.5, 42);

    SimulationEngine(make_config(bed_path1, N_SNPS, 0.5, 0.0, 123)).run();
    SimulationEngine(make_config(bed_path2, N_SNPS, 0.5, 0.0, 123)).run();

    auto phen1 = fs::path(bed_path1).replace_extension(".phen");
    auto phen2 = fs::path(bed_path2).replace_extension(".phen");
    REQUIRE(read_file_content(phen1) == read_file_content(phen2));

    auto causal1 = fs::path(bed_path1).replace_extension(".causal");
    auto causal2 = fs::path(bed_path2).replace_extension(".causal");
    REQUIRE(read_file_content(causal1) == read_file_content(causal2));
}

TEST_CASE("SimulationEngine - dominance effects emitted", "[simulate]")
{
    BedFixture fixture;
    constexpr Eigen::Index N_SNPS = 100;
    auto [bed_path, _]
        = fixture.create_bed_files(50, N_SNPS, 0.0, 0.05, 0.5, 42);

    SimulationEngine(make_config(bed_path, N_SNPS, 0.5, 0.2)).run();

    auto causal_path = fs::path(bed_path).replace_extension(".causal");
    auto effects = parse_causal_effects(causal_path, /*has_dominance=*/true);
    bool has_nonzero_dominance = std::any_of(
        effects.begin(),
        effects.end(),
        [](const auto& kv) { return std::abs(kv.second.dominance) > 1e-10; });

    REQUIRE(has_nonzero_dominance);
}

TEST_CASE("SimulationEngine - additive variance", "[simulate]")
{
    BedFixture fixture;
    constexpr Eigen::Index N_SAMPLES = 500;
    constexpr Eigen::Index N_SNPS = 200;
    constexpr double H2 = 0.5;

    auto genotypes = generate_random_genotypes(N_SAMPLES, N_SNPS, 99);
    auto [bed_path, stored_geno]
        = fixture.create_deterministic_bed_files(genotypes);

    SimulationEngine(make_config(bed_path, N_SNPS, H2)).run();

    auto effects = parse_causal_effects(
        fs::path(bed_path).replace_extension(".causal"),
        /*has_dominance=*/false);
    auto snp_to_col
        = parse_snp_indices(fs::path(bed_path).replace_extension(".bim"));
    auto [causal_geno, add_betas, _]
        = extract_causal_columns(stored_geno, effects, snp_to_col);

    process_matrix<GeneticMode::A>(
        GenotypeProcessMethod::OrthStandardize(), causal_geno);
    Eigen::VectorXd g_a = causal_geno * add_betas;

    auto phenotypes = parse_phenotypes(
        fs::path(bed_path).replace_extension(".phen"), N_SAMPLES);

    double observed_h2 = detail::var(g_a)(0) / detail::var(phenotypes)(0);
    REQUIRE_THAT(observed_h2, WithinAbs(H2, VARIANCE_TOLERANCE));
}

TEST_CASE("SimulationEngine - additive and dominance variance", "[simulate]")
{
    BedFixture fixture;
    constexpr Eigen::Index N_SAMPLES = 500;
    constexpr Eigen::Index N_SNPS = 200;
    constexpr double H2 = 0.4;
    constexpr double D2 = 0.2;

    auto genotypes = generate_random_genotypes(N_SAMPLES, N_SNPS, 99);
    auto [bed_path, stored_geno]
        = fixture.create_deterministic_bed_files(genotypes);

    SimulationEngine(make_config(bed_path, N_SNPS, H2, D2)).run();

    auto effects = parse_causal_effects(
        fs::path(bed_path).replace_extension(".causal"),
        /*has_dominance=*/true);
    auto snp_to_col
        = parse_snp_indices(fs::path(bed_path).replace_extension(".bim"));
    auto [causal_geno, add_betas, dom_betas]
        = extract_causal_columns(stored_geno, effects, snp_to_col);

    Eigen::MatrixXd x_add = causal_geno;
    process_matrix<GeneticMode::A>(
        GenotypeProcessMethod::OrthStandardize(), x_add);
    Eigen::VectorXd g_a = x_add * add_betas;

    Eigen::MatrixXd x_dom = causal_geno;
    process_matrix<GeneticMode::D>(
        GenotypeProcessMethod::OrthStandardize(), x_dom);
    Eigen::VectorXd g_d = x_dom * dom_betas;

    auto phenotypes = parse_phenotypes(
        fs::path(bed_path).replace_extension(".phen"), N_SAMPLES);

    double var_ga = detail::var(g_a)(0);
    double var_gd_raw = detail::var(g_d)(0);
    double scale = std::sqrt(var_ga * D2 / H2 / var_gd_raw);
    Eigen::VectorXd g_d_scaled = g_d * scale;

    double var_phen = detail::var(phenotypes)(0);
    REQUIRE_THAT(var_ga / var_phen, WithinAbs(H2, VARIANCE_TOLERANCE));
    REQUIRE_THAT(
        detail::var(g_d_scaled)(0) / var_phen,
        WithinAbs(D2, VARIANCE_TOLERANCE));
}
