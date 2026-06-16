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
#include <cstddef>
#include <filesystem>
#include <fstream>
#include <iterator>
#include <optional>
#include <random>
#include <sstream>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

#include <fmt/format.h>
#include <Eigen/Core>
#include <catch2/catch_test_macros.hpp>

#include "bed_fixture.h"
#include "gelex/data/genotype_method.h"
#include "gelex/data/reader.h"
#include "gelex/data/sample_id.h"
#include "gelex/io/detail/text_writer.h"
#include "gelex/simulate/effect_sampler.h"
#include "gelex/simulate/genetic_value_calculator.h"
#include "gelex/simulate/genetic_value_scaler.h"
#include "gelex/simulate/sim_types.h"
#include "gelex/types/genetic_effect_type.h"

namespace fs = std::filesystem;

using namespace gelex;  // NOLINT
using gelex::test::BedFixture;

namespace
{

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

struct SimulateScheme
{
    double heritability;
    std::vector<EffectSize> effect_sizes;
};

struct SimulateConfig
{
    int seed;
    std::string bfile_prefix;
    std::string output_prefix;
    GenotypeMethod geno_method;
    std::optional<SimulateScheme> additive;
    std::optional<SimulateScheme> dominance;
};

auto make_config(
    const fs::path& bed_path,
    Eigen::Index n_snps,
    double h2,
    double d2,
    int seed = 42) -> SimulateConfig
{
    std::optional<SimulateScheme> additive;
    if (h2 > 0.0)
    {
        additive = SimulateScheme{
            .heritability = h2,
            .effect_sizes = {{n_snps, 1.0}},
        };
    }
    std::optional<SimulateScheme> dominance;
    if (d2 > 0.0)
    {
        dominance = SimulateScheme{
            .heritability = d2,
            .effect_sizes = {{n_snps, 1.0}},
        };
    }
    return {
        .seed = seed,
        .bfile_prefix = bed_path,
        .output_prefix = bed_path,
        .geno_method = GenotypeMethod::OrthStandardize,
        .additive = std::move(additive),
        .dominance = std::move(dominance),
    };
}

auto run_simulation(const SimulateConfig& config) -> void
{
    std::mt19937_64 rng(config.seed);

    auto bim = read_bim(config.bfile_prefix + ".bim");
    auto fam = read_fam(config.bfile_prefix + ".fam");

    GeneticValueCalculator calculator(config.bfile_prefix, bim, fam);

    auto all_ids = bim.index().keys();
    std::vector<std::string_view> shuffled_ids(all_ids.begin(), all_ids.end());
    std::ranges::shuffle(shuffled_ids, rng);

    std::optional<GeneticValues> additive;
    if (config.additive)
    {
        additive.emplace();
        additive->causal_snps
            = NormalSampler(shuffled_ids, config.additive->effect_sizes)(rng);
        additive->gebv = calculator.calculate<GeneticMode::A>(
            *additive, config.geno_method);
    }

    std::optional<GeneticValues> dominance;
    if (config.dominance)
    {
        dominance.emplace();
        dominance->causal_snps
            = NormalSampler(shuffled_ids, config.dominance->effect_sizes)(rng);
        dominance->gebv = calculator.calculate<GeneticMode::D>(
            *dominance, config.geno_method);
    }

    const Eigen::Index n_samples
        = additive ? additive->gebv.size() : dominance->gebv.size();
    std::normal_distribution<double> normal01(0.0, 1.0);
    Eigen::VectorXd residual = Eigen::VectorXd::NullaryExpr(
        n_samples, [&] { return normal01(rng); });

    auto h2_of = [](const auto& scheme) -> std::optional<double>
    {
        return scheme ? std::optional<double>(scheme->heritability)
                      : std::nullopt;
    };
    GeneticValueScaler scaler(h2_of(config.additive), h2_of(config.dominance));
    scaler.scale(
        additive ? &*additive : nullptr,
        dominance ? &*dominance : nullptr,
        residual);

    Eigen::VectorXd phenotypes = additive ? additive->gebv : dominance->gebv;
    if (additive && dominance)
    {
        phenotypes += dominance->gebv;
    }
    phenotypes += residual;

    io::detail::TextWriter writer(config.output_prefix + ".phen");
    writer.write_header({"FID", "IID", "Phenotype"});
    for (Eigen::Index i = 0; i < n_samples; ++i)
    {
        auto [fid, iid] = split_sample_id(fam.index().keys()[i]);
        writer.write(fmt::format("{}\t{}\t{}", fid, iid, phenotypes(i)));
    }

    io::detail::TextWriter effect_writer(config.output_prefix + ".causal");
    auto write_single
        = [&](std::string_view column, const std::vector<CausalSnp>& snps)
    {
        effect_writer.write_header({"id", std::string(column)});
        for (const auto& snp : snps)
        {
            effect_writer.write(fmt::format("{}\t{}", snp.id, snp.effect));
        }
    };
    if (additive && dominance)
    {
        effect_writer.write_header({"id", "additive", "dominance"});
        const auto& add_snps = additive->causal_snps;
        const auto& dom_snps = dominance->causal_snps;
        for (std::size_t i = 0; i < add_snps.size(); ++i)
        {
            effect_writer.write(
                fmt::format(
                    "{}\t{}\t{}",
                    add_snps[i].id,
                    add_snps[i].effect,
                    dom_snps[i].effect));
        }
    }
    else if (additive)
    {
        write_single("additive", additive->causal_snps);
    }
    else
    {
        write_single("dominance", dominance->causal_snps);
    }
}

}  // namespace

TEST_CASE("Simulation command dataflow writes output files per mode", "[simulate]")
{
    BedFixture fixture;
    constexpr Eigen::Index N_SAMPLES = 50;
    constexpr Eigen::Index N_SNPS = 100;
    auto [bed_path, _]
        = fixture.create_bed_files(N_SAMPLES, N_SNPS, 0.0, 0.05, 0.5, 42);

    auto phen_path = fs::path(bed_path).replace_extension(".phen");
    auto causal_path = fs::path(bed_path).replace_extension(".causal");

    SECTION("A mode")
    {
        run_simulation(make_config(bed_path, N_SNPS, 0.5, 0.0));
        REQUIRE(read_first_line(phen_path) == "FID\tIID\tPhenotype");
        REQUIRE(count_lines(phen_path) == N_SAMPLES + 1);
        REQUIRE(read_first_line(causal_path) == "id\tadditive");
        REQUIRE(count_lines(causal_path) == N_SNPS + 1);
    }

    SECTION("AD mode")
    {
        run_simulation(make_config(bed_path, N_SNPS, 0.4, 0.2));
        REQUIRE(read_first_line(causal_path) == "id\tadditive\tdominance");
        REQUIRE(count_lines(causal_path) == N_SNPS + 1);
    }

    SECTION("D mode")
    {
        run_simulation(make_config(bed_path, N_SNPS, 0.0, 0.3));
        REQUIRE(read_first_line(causal_path) == "id\tdominance");
        REQUIRE(count_lines(causal_path) == N_SNPS + 1);
    }
}

TEST_CASE("Simulation command dataflow is reproducible", "[simulate]")
{
    BedFixture fixture1;
    BedFixture fixture2;
    constexpr Eigen::Index N_SNPS = 100;
    auto [bed_path1, _1]
        = fixture1.create_bed_files(50, N_SNPS, 0.0, 0.05, 0.5, 42);
    auto [bed_path2, _2]
        = fixture2.create_bed_files(50, N_SNPS, 0.0, 0.05, 0.5, 42);

    run_simulation(make_config(bed_path1, N_SNPS, 0.5, 0.0, 123));
    run_simulation(make_config(bed_path2, N_SNPS, 0.5, 0.0, 123));

    auto phen1 = fs::path(bed_path1).replace_extension(".phen");
    auto phen2 = fs::path(bed_path2).replace_extension(".phen");
    REQUIRE(read_file_content(phen1) == read_file_content(phen2));
}
