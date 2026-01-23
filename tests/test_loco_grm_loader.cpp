#include <Eigen/Core>
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include "../src/data/grm_bin_writer.h"
#include "../src/data/grm_id_writer.h"
#include "file_fixture.h"
#include "gelex/data/loco_grm_loader.h"

namespace fs = std::filesystem;
using gelex::detail::GrmBinWriter;
using gelex::detail::GrmIdWriter;
using gelex::test::FileFixture;

namespace
{
struct GrmFiles
{
    fs::path prefix;
    auto create(
        const Eigen::MatrixXd& m,
        const std::vector<std::string>& ids,
        double k) const -> void
    {
        {
            GrmBinWriter writer(prefix.string() + ".grm.bin");
            writer.write(m, k);
        }
        {
            GrmIdWriter writer(prefix.string() + ".grm.id");
            writer.write(ids);
        }
    }
};
}  // namespace

TEST_CASE("LocoGRMLoader - Basic Calculation", "[data][grm][loco]")
{
    FileFixture fixture;
    auto tmp_dir = fixture.generate_random_file_path("loco_test");
    fs::create_directories(tmp_dir);

    Eigen::Index n = 3;
    std::vector<std::string> ids = {"F1_I1", "F1_I2", "F1_I3"};

    // G_whole = (X_w * X_w') / K_w
    // G_i = (X_i * X_i') / K_i
    // G_loco = (G_w * K_w - G_i * K_i) / (K_w - K_i)

    Eigen::MatrixXd x_w = Eigen::MatrixXd::Random(n, 10);
    Eigen::MatrixXd x_i = Eigen::MatrixXd::Random(n, 3);

    double k_w = 10.0;
    double k_i = 3.0;

    Eigen::MatrixXd raw_g_w = x_w * x_w.transpose();
    Eigen::MatrixXd raw_g_i = x_i * x_i.transpose();

    GrmFiles whole_files{tmp_dir / "whole"};
    whole_files.create(raw_g_w, ids, k_w);

    GrmFiles chr_files{tmp_dir / "chr1"};
    chr_files.create(raw_g_i, ids, k_i);

    std::unordered_map<std::string, Eigen::Index> id_map;
    for (Eigen::Index i = 0; i < n; ++i)
    {
        id_map[ids[i]] = i;
    }

    gelex::LocoGRMLoader loco_loader(whole_files.prefix, id_map);
    Eigen::MatrixXd loco_grm
        = loco_loader.load_loco_grm(chr_files.prefix, id_map);

    Eigen::MatrixXd expected
        = (x_w * x_w.transpose() - x_i * x_i.transpose()) / (k_w - k_i);

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

TEST_CASE("LocoGRMLoader - Filtered Loading", "[data][grm][loco]")
{
    FileFixture fixture;
    auto tmp_dir = fixture.generate_random_file_path("loco_test_filtered");
    fs::create_directories(tmp_dir);

    Eigen::Index n = 3;
    std::vector<std::string> ids = {"F1_I1", "F1_I2", "F1_I3"};

    Eigen::MatrixXd x_w = Eigen::MatrixXd::Random(n, 10);
    Eigen::MatrixXd x_i = Eigen::MatrixXd::Random(n, 3);

    double k_w = 10.0;
    double k_i = 3.0;

    Eigen::MatrixXd raw_g_w = x_w * x_w.transpose();
    Eigen::MatrixXd raw_g_i = x_i * x_i.transpose();

    GrmFiles whole_files{tmp_dir / "whole"};
    whole_files.create(raw_g_w, ids, k_w);

    GrmFiles chr_files{tmp_dir / "chr1"};
    chr_files.create(raw_g_i, ids, k_i);

    // Test with subset of IDs and reordering
    std::unordered_map<std::string, Eigen::Index> id_map
        = {{"F1_I3", 0}, {"F1_I1", 1}};

    gelex::LocoGRMLoader loco_loader(whole_files.prefix, id_map);

    Eigen::MatrixXd loco_grm
        = loco_loader.load_loco_grm(chr_files.prefix, id_map);

    REQUIRE(loco_grm.rows() == 2);
    REQUIRE(loco_grm.cols() == 2);

    Eigen::MatrixXd x_w_subset(2, 10);
    x_w_subset.row(0) = x_w.row(2);  // I3
    x_w_subset.row(1) = x_w.row(0);  // I1

    Eigen::MatrixXd x_i_subset(2, 3);
    x_i_subset.row(0) = x_i.row(2);  // I3
    x_i_subset.row(1) = x_i.row(0);  // I1

    Eigen::MatrixXd expected = (x_w_subset * x_w_subset.transpose()
                                - x_i_subset * x_i_subset.transpose())
                               / (k_w - k_i);

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
