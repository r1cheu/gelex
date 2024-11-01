#pragma once
#include <armadillo>

namespace chenx {
using namespace arma;
template <typename T>
void value_impute(Mat<T>& genotype, const Col<T>& values);

template <typename T>
Col<T> mean_impute(Mat<T>& genotype);

template <typename T>
Col<T> median_impute(Mat<T>& genotype);
}  // namespace chenx

#include "impute_impl.h"
