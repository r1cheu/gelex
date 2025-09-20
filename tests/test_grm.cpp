#include <cstdint>

#include <Eigen/Core>
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include "gelex/data/grm.h"
//
// TEST_CASE("Genetic Relationship Matrix Test", "[Grm]")
// {
//     static constexpr double ABS_DIFF = 1e-5;
//     static constexpr uint64_t BIG_CHUNK_SIZE = 2;
//     static constexpr uint64_t SMALL_CHUNK_SIZE = 10;
//
//     using Eigen::MatrixXd;
//     using Eigen::VectorXd;
//     const std::string train_bed = std::string(GELEX_TESTS_DIR) +
//     "/data/train";
//
//     SECTION("Initialization")
//     {
//         REQUIRE_NOTHROW(gelex::GRM(train_bed, BIG_CHUNK_SIZE));
//     }
//
//     SECTION("Small chunk size")
//     {
//         gelex::GRM grm_maker(train_bed, SMALL_CHUNK_SIZE);
//         MatrixXd grm = grm_maker.compute();
//         REQUIRE(!grm.isZero());
//     }
//     SECTION("Big chunk size")
//     {
//         gelex::GRM grm_maker(train_bed, BIG_CHUNK_SIZE);
//         MatrixXd grm = grm_maker.compute();
//         REQUIRE(!grm.isZero());
//     }
//
//     SECTION("Result matches under different chunk size")
//     {
//         gelex::GRM small_maker(train_bed, SMALL_CHUNK_SIZE);
//         MatrixXd small_grm = small_maker.compute();
//
//         gelex::GRM big_maker(train_bed, BIG_CHUNK_SIZE);
//         MatrixXd big_grm = big_maker.compute();
//         REQUIRE(small_grm.isApprox(big_grm, ABS_DIFF));
//     }
//
//     SECTION("Check AddGrm Result")
//     {
//         gelex::GRM grm_maker(train_bed, BIG_CHUNK_SIZE);
//         MatrixXd grm = grm_maker.compute();
//         MatrixXd expected_grm{
//             {0.33333337, -0.33333331, 1.1589792e-08},
//             {-0.33333331, 1.5, -1.1666666},
//             {1.1589792e-08, -1.1666666, 1.1666666}};
//         REQUIRE(grm.isApprox(expected_grm, ABS_DIFF));
//     }
//
//     SECTION("Check AddGrm Center")
//     {
//         gelex::GRM grm_maker(train_bed, BIG_CHUNK_SIZE);
//         MatrixXd grm = grm_maker.compute();
//         VectorXd expected{{1.0, 0.3333333, 1.3333333, 0.6666667}};
//         REQUIRE(grm_maker.p_major().isApprox(expected / 2, ABS_DIFF));
//     }
//
//     SECTION("Check AddGrm Scale factor")
//     {
//         gelex::GRM grm_maker(train_bed, BIG_CHUNK_SIZE);
//         MatrixXd grm = grm_maker.compute();
//         REQUIRE_THAT(
//             grm_maker.scale_factor(), Catch::Matchers::WithinAbs(2,
//             ABS_DIFF));
//     }
//
//     SECTION("Check DomGrm Result")
//     {
//         gelex::GRM grm_maker(train_bed, 2);
//         MatrixXd grm = grm_maker.compute(false);
//         MatrixXd expected_grm{
//             {0.88235295, 0.35294119, -0.52941173},
//             {0.35294119, 0.88235295, -2.6490952e-08},
//             {-0.52941173, -2.6490952e-08, 1.2352941}};
//         REQUIRE(grm.isApprox(expected_grm, ABS_DIFF));
//     }
//
//     SECTION("Check DomGrm Center")
//     {
//         gelex::GRM grm_maker(train_bed, BIG_CHUNK_SIZE);
//         MatrixXd grm = grm_maker.compute(false);
//         VectorXd grm_center = grm_maker.p_major();
//         grm_center = 2 * grm_center.array() * (1 - grm_center.array());
//         VectorXd center{{0.5, 0.2777778, 0.44444442, 0.44444442}};
//         REQUIRE(grm_center.isApprox(center, ABS_DIFF));
//     }
//
//     SECTION("Check DomGrm Scale factor")
//     {
//         gelex::GRM grm_maker(train_bed, 2);
//         MatrixXd grm = grm_maker.compute(false);
//         REQUIRE_THAT(
//             grm_maker.scale_factor(),
//             Catch::Matchers::WithinAbs(0.9444445, ABS_DIFF));
//     }
// }
