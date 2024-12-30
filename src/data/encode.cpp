#include "chenx/data/encode.h"

#include <omp.h>
#include <array>

namespace chenx
{
using arma::dvec;
dmat ComputeHybirdValue(const dmat& genotype, const dvec& phenotype)
{
    dmat hybird_value_mat(2, genotype.n_cols);

#pragma omp parallel for schedule(static, 8)
    for (size_t i = 0; i < genotype.n_cols; i++)
    {
        std::array<double, 3> sum{0, 0, 0};
        std::array<size_t, 3> count{0, 0, 0};

        for (size_t j = 0; j < genotype.n_rows; j++)
        {
            if (!std::isnan(phenotype(j)))
            {
                int g = static_cast<int>(genotype(j, i));
                if (g >= 0 && g <= 2)
                {
                    sum[g] += phenotype(j);
                    count[g]++;
                }
            }
        }

        if (count[0] == 0 || count[1] == 0 || count[2] == 0)
        {
            hybird_value_mat.col(i) = {0, 1};
            continue;
        }

        std::array<double, 3> mean
            = {sum[0] / static_cast<double>(count[0]),
               sum[1] / static_cast<double>(count[1]),
               sum[2] / static_cast<double>(count[2])};
        if (mean[0] > mean[2])
        {
            double d = 2 * (mean[1] - mean[2]) / (mean[0] - mean[2]);
            hybird_value_mat.col(i) = {2, std::max(d, 0.0)};
        }
        else if (mean[0] < mean[2])
        {
            double d = 2 * (mean[1] - mean[0]) / (mean[2] - mean[0]);
            hybird_value_mat.col(i) = {0, std::max(d, 0.0)};
        }
    }
    return hybird_value_mat;
}

void HybridEncode(dmat& genotype, dmat& hybird_value)
{
#pragma omp parallel for schedule(static, 8)
    for (size_t i = 0; i < genotype.n_cols; i++)
    {
        double value = hybird_value(1, i);
        if (hybird_value(0, i) == 0)
        {  // if 0, only replace 1
            for (size_t j = 0; j < genotype.n_rows; j++)
            {
                if (genotype(j, i) == 1)
                {
                    genotype(j, i) = value;
                }
            }
        }
        else
        {  // if 2, we swap 0 and 2
            for (size_t j = 0; j < genotype.n_rows; j++)
            {
                switch (static_cast<int>(genotype(j, i)))
                {
                    case 0:
                        genotype(j, i) = 2;
                        break;
                    case 1:
                        genotype(j, i) = value;
                        break;
                    case 2:
                        genotype(j, i) = 0;
                        break;
                    default:
                        break;
                }
            }
        }
    }
}
}  // namespace chenx
