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
#include <string>
#include <utility>
#include <vector>

#include "gelex/data/dataframe/index.h"
#include "gelex/freq/design.h"

#include "cli/random_design_data.h"
#include "cli/reml_data.h"
#include "file_fixture.h"
#include "sample_id_fixture.h"

TEST_CASE(
    "Random design data serves raw and REML consumers",
    "[cli][random_design_data]")
{
    gelex::test::FileFixture files;
    const auto drand = files.create_named_text_file(
        "random_factors.tsv",
        "FID\tIID\tGroup\n"
        "F1\tI1\tA\n"
        "F1\tI2\tB\n"
        "F1\tI3\tA\n");
    const auto qrand = files.create_named_text_file(
        "random_slopes.tsv",
        "FID\tIID\tSlope\n"
        "F1\tI3\t3.0\n"
        "F1\tI1\t1.0\n"
        "F1\tI2\t2.0\n");
    const cli::RandomDesignDataConfig random_config{
        .drand_path = drand.string(), .qrand_paths = {qrand.string()}};
    const gelex::DataFrameIndex<std::string> target{std::vector<std::string>{
        gelex::make_sample_id("F1", "I2"), gelex::make_sample_id("F1", "I3")}};

    SECTION("raw loader aligns both input kinds")
    {
        cli::RandomDesignDataLoader loader(random_config);
        std::vector<const gelex::DataFrameIndex<std::string>*> indices{&target};

        loader.load_indices(indices);
        const auto names = loader.effect_names();
        REQUIRE(names.size() == 2);
        REQUIRE(names[0] == "Group");
        REQUIRE(names[1] == qrand.stem().string());
        const auto common = gelex::intersect<std::string>(indices);
        loader.gather(common);
        auto data = std::move(loader).results();

        REQUIRE(data.discrete.has_value());
        REQUIRE(data.discrete->rows() == 2);
        REQUIRE(data.discrete->col(0).as<std::string>()[0] == "B");
        REQUIRE(data.discrete->col(0).as<std::string>()[1] == "A");
        REQUIRE(data.quantitative.size() == 1);
        REQUIRE(data.quantitative[0].name == qrand.stem().string());
        REQUIRE(data.quantitative[0].frame.rows() == 2);
        REQUIRE(data.quantitative[0].frame.col(0).as<double>()[0] == 2.0);
        REQUIRE(data.quantitative[0].frame.col(0).as<double>()[1] == 3.0);
    }

    SECTION("REML loader converts shared data to kernels")
    {
        const cli::RemlDataConfig config{.designs = random_config};
        cli::RemlDataLoader loader(config);
        std::vector<const gelex::DataFrameIndex<std::string>*> indices{&target};

        loader.load_indices(indices);
        const auto common = gelex::intersect<std::string>(indices);
        loader.gather(common);
        auto designs = std::move(loader).results();

        REQUIRE(designs.size() == 2);
        REQUIRE(designs[0].name == "Group");
        REQUIRE(designs[0].kind == gelex::freq::RandomKind::Discrete);
        REQUIRE(designs[0].K.isApprox(Eigen::MatrixXd::Identity(2, 2)));
        REQUIRE(designs[1].name == qrand.stem().string());
        REQUIRE(designs[1].kind == gelex::freq::RandomKind::Quantitative);
        REQUIRE(designs[1].K.isApprox(Eigen::MatrixXd{{4.0, 6.0}, {6.0, 9.0}}));
    }
}
