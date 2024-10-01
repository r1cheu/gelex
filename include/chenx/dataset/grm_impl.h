#pragma once
#include "grm.h"

namespace chenx {
namespace dataset {
template <typename T>
void normalize(arma::Mat<T>& genotype, std::string_view method) {
    arma::Row<T> pA = arma::mean(genotype, 0) / 2;
    arma::Row<T> pa = 1 - pA;
    if (method == "add" || method == "hybrid") {
        genotype.each_row() -= 2 * pA;
    } else if (method == "dom") {
        genotype.each_row() -= 2 * (pA % pa);
    } else {
        throw std::invalid_argument("method must be 'add', 'dom' or 'hybrid'");
    }
}

template <typename T>
arma::Mat<T> cal_grm(const arma::Mat<T>& genotype) {
    arma::Mat<T> W = genotype * genotype.t();
    return W / arma::trace(W) * W.n_rows;
}

template <typename T>
arma::Mat<T> cal_grm_block(const arma::Mat<T>& genotype, const arma::uword block_size) {
    auto n = genotype.n_rows;
    arma::Mat<T> W(n, n, arma::fill::zeros);

    for (size_t i = 0; i < genotype.n_cols; i += block_size) {
        auto colEnd = std::min(i + block_size, genotype.n_cols);
        arma::Mat<T> genotypeBlock = genotype.cols(i, colEnd - 1);
        W += genotypeBlock * genotypeBlock.t();
    }
    W = W / arma::trace(W) * W.n_rows;

    return W;
}
}  // namespace dataset
}  // namespace chenx
