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
#include <ranges>
#include <string>
#include <utility>
#include <vector>

#include "gelex/data/bed.h"
#include "gelex/data/dataframe/index.h"

#include "bed_fixture.h"
#include "cli/mcmc/data.h"
#include "cli/random_design_data.h"
#include "sample_id_fixture.h"

TEST_CASE(
    "MCMC data loader aligns Bed and builds explicit random designs",
    "[cli][mcmc][data]")
{
    gelex::test::BedFixture fixture;
    const auto [prefix, genotypes] = fixture.create_deterministic_bed_files(
        Eigen::MatrixXd{{0.0}, {1.0}, {2.0}}, {"I1", "I2", "I3"});
    static_cast<void>(genotypes);
    auto& files = fixture.get_file_fixture();
    const auto drand = files.create_named_text_file(
        "random_factors.tsv",
        "FID\tIID\tGroup\n"
        "fam3\tI3\tA\n"
        "fam1\tI1\tA\n"
        "fam2\tI2\tB\n");
    const auto qrand = files.create_named_text_file(
        "random_slopes.tsv",
        "FID\tIID\tSlope\n"
        "fam1\tI1\t1.0\n"
        "fam3\tI3\t3.0\n"
        "fam2\tI2\t2.0\n");
    const cli::RandomDesignDataConfig config{
        .drand_path = drand.string(), .qrand_paths = {qrand.string()}};
    cli::McmcDataLoader loader(gelex::open_bed(prefix.string()), config);
    const gelex::DataFrameIndex<std::string> target{std::vector<std::string>{
        gelex::make_sample_id("fam2", "I2"),
        gelex::make_sample_id("fam3", "I3")}};
    std::vector<const gelex::DataFrameIndex<std::string>*> indices{&target};

    loader.load_indices(indices);
    const auto common = gelex::intersect<std::string>(indices);
    loader.gather(common);
    auto data = std::move(loader).results();

    REQUIRE(std::ranges::equal(data.bed.sample_index().keys(), target.keys()));
    REQUIRE(data.random.size() == 2);
    REQUIRE(data.random[0].name() == "Group");
    REQUIRE(
        data.random[0].X().isApprox(Eigen::MatrixXd{{0.0, 1.0}, {1.0, 0.0}}));
    REQUIRE(data.random[1].name() == qrand.stem().string());
    REQUIRE(data.random[1].X().isApprox(Eigen::MatrixXd{{2.0}, {3.0}}));
}
