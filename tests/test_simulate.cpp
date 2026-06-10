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
#include <filesystem>
#include <fstream>
#include <iterator>
#include <optional>
#include <sstream>
#include <string>
#include <utility>

#include <Eigen/Core>
#include <catch2/catch_test_macros.hpp>

#include "bed_fixture.h"
#include "gelex/data/genotype/method.h"
#include "gelex/engine/simulation.h"

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

auto make_config(
    const fs::path& bed_path,
    Eigen::Index n_snps,
    double h2,
    double d2,
    int seed = 42) -> SimulationEngine::Config
{
    std::optional<SimulationEngine::SimulateScheme> additive;
    if (h2 > 0.0)
    {
        additive = SimulationEngine::SimulateScheme{
            .heritability = h2,
            .effect_sizes = {{n_snps, 1.0}},
        };
    }
    std::optional<SimulationEngine::SimulateScheme> dominance;
    if (d2 > 0.0)
    {
        dominance = SimulationEngine::SimulateScheme{
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

}  // namespace

TEST_CASE("SimulationEngine - output files per mode", "[simulate]")
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
        SimulationEngine(make_config(bed_path, N_SNPS, 0.5, 0.0)).run();
        REQUIRE(read_first_line(phen_path) == "FID\tIID\tPhenotype");
        REQUIRE(count_lines(phen_path) == N_SAMPLES + 1);
        REQUIRE(read_first_line(causal_path) == "id\tadditive");
        REQUIRE(count_lines(causal_path) == N_SNPS + 1);
    }

    SECTION("AD mode")
    {
        SimulationEngine(make_config(bed_path, N_SNPS, 0.4, 0.2)).run();
        REQUIRE(read_first_line(causal_path) == "id\tadditive\tdominance");
        REQUIRE(count_lines(causal_path) == N_SNPS + 1);
    }

    SECTION("D mode")
    {
        SimulationEngine(make_config(bed_path, N_SNPS, 0.0, 0.3)).run();
        REQUIRE(read_first_line(causal_path) == "id\tdominance");
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
}
