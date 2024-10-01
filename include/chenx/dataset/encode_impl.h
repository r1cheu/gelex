#pragma once
#include <omp.h>
#include "encode.h"

namespace chenx {
namespace dataset {
template <typename T>
void dominance(arma::Mat<T>& genotype) {
    genotype.replace(T(2), T(0));
}

template <typename T>
arma::Mat<T> hybird_value(const arma::Mat<T>& genotype, const arma::Col<T>& phenotype) {
    arma::Mat<T> hybird_value_mat(2, genotype.n_cols);
#pragma omp parallel for
    for (size_t i = 0; i < genotype.n_cols; i++) {
        double sum0 = 0, sum1 = 0, sum2 = 0;
        size_t elems_count0 = 0, elems_count1 = 0, elems_count2 = 0;

        for (size_t j = 0; j < genotype.n_rows; j++) {
            if (genotype(j, i) == 0 && !std::isnan(phenotype(j))) {
                sum0 += phenotype(j);
                elems_count0++;
            } else if (genotype(j, i) == 1 && !std::isnan(phenotype(j))) {
                sum1 += phenotype(j);
                elems_count1++;
            } else if (genotype(j, i) == 2 && !std::isnan(phenotype(j))) {
                sum2 += phenotype(j);
                elems_count2++;
            }
        }

        if (elems_count0 == 0 || elems_count1 == 0 || elems_count2 == 0) {
            // TODO add user warning
            hybird_value_mat.col(i) = {0, 1};
            continue;
        }

        double mean0 = sum0 / elems_count0;
        double mean1 = sum1 / elems_count1;
        double mean2 = sum2 / elems_count2;

        if (mean0 > mean2) {
            double d = 2 * (mean1 - mean2) / (mean0 - mean2);
            if (d < 0) {
                d = 0;
            }
            hybird_value_mat.col(i) = {2, static_cast<T>(d)};
        }

        if (mean0 < mean2) {
            double d = 2 * (mean1 - mean0) / (mean2 - mean0);
            if (d < 0) {
                d = 0;
            }
            hybird_value_mat.col(i) = {0, static_cast<T>(d)};
        }
    }
    return hybird_value_mat;
}

template <typename T>
void hybird(arma::Mat<T>& genotype, const arma::Mat<T>& hybird_value) {
#pragma omp parallel for
    for (size_t i = 0; i < genotype.n_cols; i++) {
        T value = hybird_value(1, i);
        if (hybird_value(0, i) == 0) {  // if 0, only replace 1
            for (size_t j = 0; j < genotype.n_rows; j++) {
                if (genotype(j, i) == 1) {
                    genotype(j, i) = value;
                }
            }
        } else {  // if 2, we swap 0 and 2
            for (size_t j = 0; j < genotype.n_rows; j++) {
                if (genotype(j, i) == 0) {
                    genotype(j, i) = 2;
                } else if (genotype(j, i) == 1) {
                    genotype(j, i) = value;
                } else if (genotype(j, i) == 2) {
                    genotype(j, i) = 0;
                }
            }
        }
    }
}
}  // namespace dataset
}  // namespace chenx
