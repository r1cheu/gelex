#include <gtest/gtest.h>
#include "armadillo"
#include "chenx/dataset/encode.h"
#include "chenx/dataset/grm.h"
#include "chenx/dataset/impute.h"
#include "chenx/optim/zkztr.h"

namespace chenx {
using namespace arma;
TEST(FillNaTest, Mean) {
    dmat X = {{NAN, 2.0, 3.0}, {4.0, NAN, 6.0}, {7.0, 8.0, NAN}};
    mean_impute(X);
    dmat expected = {{5.5, 2.0, 3.0}, {4.0, 5.0, 6.0}, {7.0, 8.0, 4.5}};
    ASSERT_TRUE(arma::approx_equal(X, expected, "absdiff", 1e-10));
}

TEST(FillNaTest, MedianEven) {
    dmat X = {{NAN, 2.0, 3.0}, {4.0, NAN, 6.0}, {7.0, 8.0, NAN}};
    median_impute(X);

    dmat expected = {{5.5, 2.0, 3.0}, {4.0, 5.0, 6.0}, {7.0, 8.0, 4.5}};

    ASSERT_TRUE(arma::approx_equal(X, expected, "absdiff", 1e-10));
}

TEST(FillNaTest, MedianOdd) {
    dmat X = {{NAN, 2.0, 3.0}, {4.0, NAN, 6.0}, {7.0, 8.0, NAN}, {1.0, 2.0, 3.0}};
    median_impute(X);

    dmat expected = {{4.0, 2.0, 3.0}, {4.0, 2.0, 6.0}, {7.0, 8.0, 3.0}, {1.0, 2.0, 3.0}};

    ASSERT_TRUE(arma::approx_equal(X, expected, "absdiff", 1e-10));
}

TEST(EncodeTest, Hybird) {
    dmat X = {{1, 0, 2, 2}, {1, 2, 2, 1}, {2, 2, 2, 2}, {2, 2, 2, 1}, {1, 0, 2, 2}};
    dmat guide = {{0, 0, 0, 2}, {1, 1.5, 2, 2.5}};
    hybird(X, guide);

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
    dmat result = hybird_value(X, phenotype);
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
    dmat result = hybird_value(X, phenotype);
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
    dmat result = hybird_value(X, phenotype);
    dmat expected = {{0, 0, 2}, {2, 0, 3}};
    ASSERT_TRUE(approx_equal(result, expected, "absdiff", 1e-10));
}

// Helper function to create identity matrix
template <typename eT>
Mat<eT> create_identity_matrix(size_t size) {
    return eye(size, size);
}

template <typename eT>
SpMat<eT> create_identity_sparse_matrix(size_t size) {
    return speye(size, size);
}

// Helper function to create non-identity matrix
template <typename eT>
Mat<eT> create_non_identity_matrix(size_t size) {
    Mat<eT> mat(size, size, fill::randu);
    return mat;
}

// Test when both z and k are identity matrices
TEST(CalZkzTest, BothIdentity) {
    SpMat<double> z = create_identity_sparse_matrix<double>(3);
    Mat<double> k = create_identity_matrix<double>(3);
    Mat<double> result = cal_zkz(z, k);
    ASSERT_TRUE(approx_equal(result, k, "absdiff", 0.0001));
}

// Test when z is identity and k is non-identity
TEST(CalZkzTest, ZIdentityKNonIdentity) {
    SpMat<double> z = create_identity_sparse_matrix<double>(3);
    Mat<double> k = create_non_identity_matrix<double>(3);
    Mat<double> result = cal_zkz(z, k);
    ASSERT_TRUE(approx_equal(result, k, "absdiff", 0.0001));
}

// Test when z is non-identity and k is identity
TEST(CalZkzTest, ZNonIdentityKIdentity) {
    SpMat<double> z = create_identity_sparse_matrix<double>(3);
    Mat<double> k = create_identity_matrix<double>(3);
    Mat<double> result = cal_zkz(z, k);
    Mat<double> expected(z * z.t());
    ASSERT_TRUE(approx_equal(result, expected, "absdiff", 0.0001));
}

// Test when both z and k are non-identity matrices
TEST(CalZkzTest, BothNonIdentity) {
    SpMat<double> z = create_identity_sparse_matrix<double>(3);
    Mat<double> k = create_non_identity_matrix<double>(3);
    Mat<double> result = cal_zkz(z, k);
    Mat<double> expected = z * k * z.t();
    ASSERT_TRUE(approx_equal(result, expected, "absdiff", 0.0001));
}

}  // namespace chenx
