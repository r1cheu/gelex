#pragma once
#include <armadillo>

namespace chenx {
namespace dataset {
template <typename T>
void normalize(arma::Mat<T>& genotype, std::string_view method);
template <typename T>
arma::Mat<T> cal_grm(const arma::Mat<T>& genotype);
template <typename T>
arma::Mat<T> cal_grm_block(const arma::Mat<T>& genotype, const arma::uword block_size);
}  // namespace dataset
}  // namespace chenx

#include "grm_impl.h"
