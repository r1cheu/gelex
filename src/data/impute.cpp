#include <chenx/data/impute.h>

#include <cmath>
#include <cstddef>

#include <iostream>
#include <vector>
#include "omp.h"

#include <armadillo>

namespace chenx
{
dvec MeanImpute(dmat& genotype)
{
    dvec means(genotype.n_cols);

#pragma omp parallel for default(none) shared(std::cerr, genotype, means)
    for (size_t i = 0; i < genotype.n_cols; i++)
    {
        double sum = 0;
        size_t elem_count = 0;
        std::vector<size_t> index;

        for (size_t j = 0; j < genotype.n_rows; j++)
        {
            if (std::isnan(genotype(j, i)))
            {
                index.emplace_back(j);
            }
            else
            {
                sum += genotype(j, i);
                elem_count++;
            }
        }

        if (elem_count == 0)
        {
#pragma omp critical
            {
                std::cerr << "Warning: All elements are missing in column " << i
                          << ". Skipping imputation." << '\n';
            }
            continue;  // 跳过当前列的处理
        }

        double mean = sum / static_cast<double>(elem_count);
        means(i) = mean;
        for (const auto& j : index)
        {
            genotype(j, i) = mean;
        }
    }

    return means;
}

dvec MedianImpute(dmat& genotype)
{
    dvec medians(genotype.n_cols);

#pragma omp parallel for default(none) shared(std::cerr, genotype, medians)
    for (size_t i = 0; i < genotype.n_cols; i++)
    {
        std::vector<size_t> index;
        std::vector<double> elems_to_keep;

        for (size_t j = 0; j < genotype.n_rows; j++)
        {
            if (std::isnan(genotype(j, i)))
            {
                index.emplace_back(j);
            }
            else
            {
                elems_to_keep.push_back(genotype(j, i));
            }
        }

        if (elems_to_keep.size() == 0)
        {
#pragma omp critical
            {
                std::cerr << "Warning: All elements are missing in column " << i
                          << ". Skipping imputation." << '\n';
            }
            continue;  // 跳过当前列的处理
        }

        double median = arma::median(arma::vec(elems_to_keep));

        medians(i) = median;
        for (const auto& j : index)
        {
            genotype(j, i) = median;
        }
    }

    return medians;
}

void ValueImpute(dmat& genotype, const dvec& values)
{
#pragma omp parallel for default(none) shared(std::cerr, genotype, values)
    for (size_t i = 0; i < genotype.n_cols; i++)
    {
        size_t elem_count = 0;
        for (size_t j = 0; j < genotype.n_rows; j++)
        {
            if (std::isnan(genotype(j, i)))
            {
                genotype(j, i) = values(i);
            }
            else
            {
                elem_count++;
            }
        }

        if (elem_count == 0)
        {
#pragma omp critical
            {
                std::cerr << "Warning: All elements are missing in column " << i
                          << ". Skipping imputation." << '\n';
            }
            continue;  // 跳过当前列的处理
        }
    }
}
}  // namespace chenx
