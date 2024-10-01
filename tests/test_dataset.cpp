#include <gtest/gtest.h>
#include "chenx/dataset/encode.h"
#include "chenx/dataset/grm.h"
#include "chenx/dataset/impute.h"
using namespace arma;

TEST(FillNaTest, Mean) {
    dmat X = {{NAN, 2.0, 3.0}, {4.0, NAN, 6.0}, {7.0, 8.0, NAN}};
    chenx::dataset::mean_impute(X);
    dmat expected = {{5.5, 2.0, 3.0}, {4.0, 5.0, 6.0}, {7.0, 8.0, 4.5}};
    ASSERT_TRUE(arma::approx_equal(X, expected, "absdiff", 1e-10));
}

TEST(FillNaTest, MedianEven) {
    dmat X = {{NAN, 2.0, 3.0}, {4.0, NAN, 6.0}, {7.0, 8.0, NAN}};
    chenx::dataset::median_impute(X);

    dmat expected = {{5.5, 2.0, 3.0}, {4.0, 5.0, 6.0}, {7.0, 8.0, 4.5}};

    ASSERT_TRUE(arma::approx_equal(X, expected, "absdiff", 1e-10));
}

TEST(FillNaTest, MedianOdd) {
    dmat X = {{NAN, 2.0, 3.0}, {4.0, NAN, 6.0}, {7.0, 8.0, NAN}, {1.0, 2.0, 3.0}};
    chenx::dataset::median_impute(X);

    dmat expected = {{4.0, 2.0, 3.0}, {4.0, 2.0, 6.0}, {7.0, 8.0, 3.0}, {1.0, 2.0, 3.0}};

    ASSERT_TRUE(arma::approx_equal(X, expected, "absdiff", 1e-10));
}

TEST(NormalizeTest, Invalid) {
    dmat X(3, 3, fill::randu);
    ASSERT_THROW(chenx::dataset::normalize(X, "invalid_method"), std::invalid_argument);
}

TEST(NormalizeTest, NormlizeAdd) {
    dmat X = {{1, 0, 2, 2}, {1, 2, 2, 1}, {2, 2, 2, 2}, {2, 2, 2, 1}};
    chenx::dataset::normalize(X, "add");
    dmat Y = {{-0.5, -1.5, 0, 0.5}, {-0.5, 0.5, 0, -0.5}, {0.5, 0.5, 0, 0.5}, {0.5, 0.5, 0, -0.5}};
    bool result = arma::approx_equal(X, Y, "absdiff", 1e-10);
    ASSERT_TRUE(result);
}

TEST(NormalizeTest, NormlizeDom) {
    dmat X = {{1, 0, 2, 2}, {1, 2, 2, 1}, {2, 2, 2, 2}, {2, 2, 2, 1}};
    chenx::dataset::normalize(X, "dom");
    dmat Y = {{0.625, -0.375, 2, 1.625}, {0.625, 1.625, 2, 0.625}, {1.625, 1.625, 2, 1.625}, {1.625, 1.625, 2, 0.625}};
    bool result = arma::approx_equal(X, Y, "absdiff", 1e-10);
    ASSERT_TRUE(result);
}

TEST(EncodeTest, Hybird) {
    dmat X = {{1, 0, 2, 2}, {1, 2, 2, 1}, {2, 2, 2, 2}, {2, 2, 2, 1}, {1, 0, 2, 2}};
    dmat guide = {{0, 0, 0, 2}, {1, 1.5, 2, 2.5}};
    chenx::dataset::hybird(X, guide);

    dmat Y = {{1, 0, 2, 0}, {1, 2, 2, 2.5}, {2, 2, 2, 0}, {2, 2, 2, 2.5}, {1, 0, 2, 0}};
    bool result = approx_equal(X, Y, "absdiff", 1e-10);
}

TEST(HybirdValue, Basic) {
    dmat X = {
        {0, 1, 2},
        {1, 0, 2},
        {2, 1, 0},
        {1, 2, 1},
    };
    dvec phenotype = {1.0, 2.0, 3.0, 4.0};
    dmat result = chenx::dataset::hybird_value(X, phenotype);
    dmat expected = {{0, 0, 2}, {2, 0, 10 / 3.0}};
    ASSERT_TRUE(approx_equal(result, expected, "absdiff", 1e-10));
}

TEST(HybirdValue, MissGenotypeHandling) {
    dmat X = {
        {0, 1, 2},
        {1, 0, 2},
        {2, 1, 0},
        {1, 2, 0},
    };
    dvec phenotype = {1.0, 2.0, 3.0, 4.0};
    dmat result = chenx::dataset::hybird_value(X, phenotype);
    result.print();
    std::cout << std::endl;
    dmat expected = {{0, 0, 0}, {2, 0, 1}};
    expected.print();
    ASSERT_TRUE(approx_equal(result, expected, "absdiff", 1e-10));
}

TEST(HybirdValue, NaNHandling) {
    dmat X = {
        {0, 1, 2},
        {1, 0, datum::nan},
        {2, 1, 0},
        {1, 2, 1},
    };
    dvec phenotype = {1.0, 2.0, 3.0, 4.0};
    dmat result = chenx::dataset::hybird_value(X, phenotype);
    dmat expected = {{0, 0, 2}, {2, 0, 3}};
    ASSERT_TRUE(approx_equal(result, expected, "absdiff", 1e-10));
}
