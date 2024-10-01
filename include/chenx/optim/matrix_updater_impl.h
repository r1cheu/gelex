#pragma once
#include "matrix_updater.h"

namespace chenx {
namespace optim {
using namespace arma;
template <typename eT>
MatrixUpdater<eT>::MatrixUpdater(const Mat<eT>& X, const Col<eT>& y, const Cube<eT>& zkztr)
    : _X{X}, _y{y}, _zkztr{zkztr} {
    auto n_fixed = _X.n_cols;
    auto N = _X.n_rows;

    _proj_y.zeros(N);

    _v.zeros(N, N);
    _vi.zeros(N, N);
    _proj.zeros(N, N);
    _txvx.zeros(n_fixed, n_fixed);

    _pdv.zeros(N, N, zkztr.n_slices);
}

template <typename eT>
void MatrixUpdater<eT>::_cal_v(const Col<eT>& var) {
    _v = _zkztr.slice(0) * var(0);
    for (size_t i{1}; i < var.n_elem; ++i) {
        _v += _zkztr.slice(i) * var(i);
    }
}

template <typename eT>
void MatrixUpdater<eT>::_cal_proj_matrix() {
    _vi = inv_sympd(_v);
    Mat<eT> vx = _vi * _X;
    _txvx = _X.t() * vx;
    _proj = _vi - vx * solve(_txvx, vx.t());
}

template <typename eT>
void MatrixUpdater<eT>::_cal_pdv() {
    for (size_t i{0}; i < _zkztr.n_slices; ++i) {
        _pdv.slice(i) = _proj * _zkztr.slice(i);
    }
}

template <typename eT>
void MatrixUpdater<eT>::update(const Col<eT>& var) {
    _cal_v(var);
    _cal_proj_matrix();
    _cal_pdv();
    _proj_y = _proj * _y;
}
}  // namespace optim
}  // namespace chenx
