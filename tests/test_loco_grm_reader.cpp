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
#include <catch2/matchers/catch_matchers.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <filesystem>
#include <string>
#include <string_view>
#include <vector>

#include "gelex/data/dataframe/index.h"
#include "gelex/data/reader.h"
#include "gelex/data/writer.h"
#include "gelex/io/loco_reader.h"

#include "file_fixture.h"
#include "sample_id_fixture.h"

namespace fs = std::filesystem;
using gelex::test::FileFixture;
namespace df = gelex;

namespace
{
auto sid(std::string_view fid, std::string_view iid) -> std::string
{
    return gelex::make_sample_id(fid, iid);
}

struct GrmFiles
{
    fs::path prefix;
    auto create(const Eigen::MatrixXd& m, const std::vector<std::string>& ids)
        const -> void
    {
        gelex::write_grm(prefix.string(), m, ids);
    }
};
}  // namespace

TEST_CASE("LocoReader - Basic Calculation", "[data][grm][loco]")
{
    FileFixture fixture;
    auto tmp_dir = fixture.generate_random_file_path("loco_test");
    fs::create_directories(tmp_dir);

    Eigen::Index n = 3;
    std::vector<std::string> ids
        = {sid("F1", "I1"), sid("F1", "I2"), sid("F1", "I3")};

    // G_loco = (G_w - G_i) / (k_w - k_i)
    // where k = trace(G) / n

    Eigen::MatrixXd x_w = Eigen::MatrixXd::Random(n, 10);
    Eigen::MatrixXd x_i = Eigen::MatrixXd::Random(n, 3);

    Eigen::MatrixXd raw_g_w = x_w * x_w.transpose();
    Eigen::MatrixXd raw_g_i = x_i * x_i.transpose();

    // k values are now computed from trace/n
    double k_w = raw_g_w.trace() / static_cast<double>(n);
    double k_i = raw_g_i.trace() / static_cast<double>(n);

    GrmFiles whole_files{tmp_dir / "whole"};
    whole_files.create(raw_g_w, ids);

    GrmFiles chr_files{tmp_dir / "chr1"};
    chr_files.create(raw_g_i, ids);

    df::DataFrameIndex<std::string> sample_index(ids);

    Eigen::MatrixXd g_whole
        = gelex::read_grm(whole_files.prefix.string(), &sample_index, false);
    gelex::LocoReader loco_reader(g_whole);
    Eigen::MatrixXd loco_grm;
    loco_reader.load_into(chr_files.prefix, sample_index, loco_grm);

    Eigen::MatrixXd expected = (raw_g_w - raw_g_i) / (k_w - k_i);

    for (Eigen::Index i = 0; i < n; ++i)
    {
        for (Eigen::Index j = 0; j < n; ++j)
        {
            REQUIRE_THAT(
                loco_grm(i, j),
                Catch::Matchers::WithinRel(expected(i, j), 1e-5));
        }
    }
}

TEST_CASE("LocoReader - Filtered Loading", "[data][grm][loco]")
{
    FileFixture fixture;
    auto tmp_dir = fixture.generate_random_file_path("loco_test_filtered");
    fs::create_directories(tmp_dir);

    Eigen::Index n = 3;
    std::vector<std::string> ids
        = {sid("F1", "I1"), sid("F1", "I2"), sid("F1", "I3")};

    Eigen::MatrixXd x_w = Eigen::MatrixXd::Random(n, 10);
    Eigen::MatrixXd x_i = Eigen::MatrixXd::Random(n, 3);

    Eigen::MatrixXd raw_g_w = x_w * x_w.transpose();
    Eigen::MatrixXd raw_g_i = x_i * x_i.transpose();

    GrmFiles whole_files{tmp_dir / "whole"};
    whole_files.create(raw_g_w, ids);

    GrmFiles chr_files{tmp_dir / "chr1"};
    chr_files.create(raw_g_i, ids);

    // Test with subset of IDs and reordering
    df::DataFrameIndex<std::string> sample_index(
        std::vector<std::string>{sid("F1", "I3"), sid("F1", "I1")});

    Eigen::MatrixXd g_whole
        = gelex::read_grm(whole_files.prefix.string(), &sample_index, false);
    gelex::LocoReader loco_reader(g_whole);

    Eigen::MatrixXd loco_grm;
    loco_reader.load_into(chr_files.prefix, sample_index, loco_grm);

    REQUIRE(loco_grm.rows() == 2);
    REQUIRE(loco_grm.cols() == 2);

    // Compute expected values using filtered matrices
    Eigen::MatrixXd x_w_subset(2, 10);
    x_w_subset.row(0) = x_w.row(2);  // I3
    x_w_subset.row(1) = x_w.row(0);  // I1

    Eigen::MatrixXd x_i_subset(2, 3);
    x_i_subset.row(0) = x_i.row(2);  // I3
    x_i_subset.row(1) = x_i.row(0);  // I1

    Eigen::MatrixXd g_w_subset = x_w_subset * x_w_subset.transpose();
    Eigen::MatrixXd g_i_subset = x_i_subset * x_i_subset.transpose();

    // k values are computed from filtered matrices' trace/n
    double k_w = g_w_subset.trace() / 2.0;
    double k_i = g_i_subset.trace() / 2.0;

    Eigen::MatrixXd expected = (g_w_subset - g_i_subset) / (k_w - k_i);

    for (Eigen::Index i = 0; i < 2; ++i)
    {
        for (Eigen::Index j = 0; j < 2; ++j)
        {
            REQUIRE_THAT(
                loco_grm(i, j),
                Catch::Matchers::WithinRel(expected(i, j), 1e-5));
        }
    }
}
