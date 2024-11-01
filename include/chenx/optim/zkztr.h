#pragma once
#include <cmath>
#include <vector>
#include "armadillo"

namespace chenx {
using namespace arma;
template <typename eT>
bool check_identity(const Mat<eT>& inputs) {
    if (!inputs.is_square()) {
        return false;
    }

    if (!inputs.is_diagmat()) {
        return false;
    }

    for (size_t i = 0; i < inputs.n_rows; ++i) {
        if (inputs(i, i) != 1) {
            return false;
        }
    }
    return true;
}

template <typename eT>
bool check_identity(const SpMat<eT>& inputs) {
    if (!inputs.is_square()) {
        return false;
    }

    if (!inputs.is_diagmat()) {
        return false;
    }

    for (size_t i = 0; i < inputs.n_rows; ++i) {
        if (inputs(i, i) != 1) {
            return false;
        }
    }
    return true;
}

template <typename eT>
std::vector<SpMat<eT>> create_z(const uword& n_z, const Col<uword>& z_index, const uword& n) {
    std::vector<SpMat<eT>> z;
    for (size_t i = 0; i < n_z; ++i) {
        z.push_back(speye<SpMat<eT>>(n, n).cols(z_index).t());
    }
    return z;
}

template <typename eT>
std::vector<SpMat<eT>> create_z(const uword& n_z, const uword& n) {
    std::vector<SpMat<eT>> z;
    for (size_t i = 0; i < n_z; ++i) {
        z.push_back(speye<SpMat<eT>>(n, n));
    }
    return z;
}

template <typename eT>
Mat<eT> cal_zkz(const SpMat<eT>& z, const Mat<eT>& k) {
    bool z_identity = check_identity(z);
    bool k_identity = check_identity(k);

    if (k_identity) {
        return z_identity ? k : Mat<eT>(z * z.t());
    } else {
        return z_identity ? k : Mat<eT>(z * k * z.t());
    }
}

template <typename eT>
Cube<eT> cal_zkztr(const std::vector<SpMat<eT>>& z, const Cube<eT>& k) {
    auto n = z[0].n_rows;
    Cube<eT> result(n, n, k.n_slices + 1, fill::zeros);

    for (size_t i = 0; i < k.n_slices; ++i) {
        result.slice(i) = cal_zkz(z[i], k.slice(i));
    }

    result.slice(k.n_slices).eye();
    return result;
}

}  // namespace chenx
