#pragma once
#include "impute.h"
#include "omp.h"

namespace chenx {
namespace dataset {

template <typename T>
void mean_impute(arma::Mat<T>& genotype) {
#pragma omp parallel for
    for (size_t i = 0; i < genotype.n_cols; i++) {
        double sum = 0;
        size_t elem_count = 0;

        std::vector<size_t> index;
        for (size_t j = 0; j < genotype.n_rows; j++) {
            if (std::isnan(genotype(j, i))) {
                index.emplace_back(j);
            } else {
                sum += genotype(j, i);
                elem_count++;
            }
        }

        if (elem_count == 0) {
            throw std::runtime_error("All elements are missing in column " + std::to_string(i));
        }

        double mean = sum / elem_count;
        for (const auto& j : index) {
            genotype(j, i) = mean;
        }
    }
}

template <typename T>
void median_impute(arma::Mat<T>& genotype) {
#pragma omp parallel for
    for (size_t i = 0; i < genotype.n_cols; i++) {
        std::vector<size_t> index;
        std::vector<double> elems_to_keep;

        for (size_t j = 0; j < genotype.n_rows; j++) {
            if (std::isnan(genotype(j, i))) {
                index.emplace_back(j);
            } else {
                elems_to_keep.push_back(genotype(j, i));
            }
        }

        if (elems_to_keep.size() == 0) {
            throw std::runtime_error("All elements are missing in column " + std::to_string(i));
        }

        double median = arma::median(arma::vec(elems_to_keep));

        for (const auto& j : index) {
            genotype(j, i) = median;
        }
    }
}
}  // namespace dataset
}  // namespace chenx
