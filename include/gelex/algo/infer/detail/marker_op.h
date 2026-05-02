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

#ifndef GELEX_ALGO_INFER_DETAIL_MARKER_OP_H_
#define GELEX_ALGO_INFER_DETAIL_MARKER_OP_H_

#include <cassert>
#include <cstdint>
#include <span>

#ifdef USE_MKL
#include <mkl.h>
#else
#include <cblas.h>
#endif

#include <Eigen/Core>

namespace gelex::detail
{

template <typename DerivedX, typename DerivedY>
inline double blas_ddot(
    const Eigen::DenseBase<DerivedX>& x,
    const Eigen::DenseBase<DerivedY>& y)
{
    EIGEN_STATIC_ASSERT_VECTOR_ONLY(DerivedX);
    EIGEN_STATIC_ASSERT_VECTOR_ONLY(DerivedY);
    assert(x.size() == y.size() && "blas_ddot: vector sizes do not match.");

    const int n = static_cast<int>(x.size());
    const int incx = 1;
    const int incy = 1;

    return cblas_ddot(n, x.derived().data(), incx, y.derived().data(), incy);
}

template <typename DerivedX, typename DerivedY>
inline void blas_daxpy(
    double alpha,
    const Eigen::DenseBase<DerivedX>& x,
    Eigen::DenseBase<DerivedY>& y)
{
    EIGEN_STATIC_ASSERT_VECTOR_ONLY(DerivedX);
    EIGEN_STATIC_ASSERT_VECTOR_ONLY(DerivedY);
    assert(x.size() == y.size() && "blas_daxpy: vector sizes do not match.");

    const int n = static_cast<int>(x.size());
    const double a = alpha;
    const int incx = 1;
    const int incy = 1;

    cblas_daxpy(n, a, x.derived().data(), incx, y.derived().data(), incy);
}

// Describes a single marker's coefficient transition.
// old_class / new_class default to -1 (no mixture context); component classes
// are 1-based indices into component_u (class 0 is the spike at zero).
struct MarkerTransition
{
    double old_value{0.0};
    double new_value{0.0};
    int8_t old_class{-1};
    int8_t new_class{-1};
};

// Multi-version dispatched (default / AVX2 / AVX-512F) inner kernel for the
// fused marker update; defined in marker_op.cpp.
auto apply_marker_update_impl(
    double* y_adj,
    double* u,
    std::span<Eigen::VectorXd> component_u,
    const double* col,
    Eigen::Index n,
    const MarkerTransition& tx) -> void;

// Fused single-pass marker update:
//   y_adj += diff * col
//   u     -= diff * col           (diff = old_value - new_value)
//   component_u[old_class-1] -= old_value * col   (if old_class > 0)
//   component_u[new_class-1] += new_value * col   (if new_class > 0)
// Pass an empty span when no mixture tracking is needed.
template <typename DerivedCol>
inline auto apply_marker_update(
    Eigen::Ref<Eigen::VectorXd> y_adj,
    Eigen::Ref<Eigen::VectorXd> u,
    std::span<Eigen::VectorXd> component_u,
    const Eigen::DenseBase<DerivedCol>& col,
    const MarkerTransition& tx) -> void
{
    EIGEN_STATIC_ASSERT_VECTOR_ONLY(DerivedCol);
    assert(col.size() == y_adj.size());
    assert(col.size() == u.size());
    apply_marker_update_impl(
        y_adj.data(),
        u.data(),
        component_u,
        col.derived().data(),
        col.size(),
        tx);
}

}  // namespace gelex::detail

#endif  // GELEX_ALGO_INFER_DETAIL_MARKER_OP_H_
