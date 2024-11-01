#pragma once
#include "omp.h"

namespace chenx
{
template <typename KernelType, typename MatType>
class NaiveKernelRule
{
   public:
    /**
     * Construct the exact kernel matrix.
     *
     * @param data Input data points.
     */
    static MatType ApplyKernelMatrix(const MatType& data, KernelType kernel = KernelType())
    {
        // Construct the kernel matrix.
        MatType kernelMatrix;
        // Resize the kernel matrix to the right size.
        kernelMatrix.set_size(data.n_cols, data.n_cols);

        // Note that we only need to calculate the upper triangular part of the
        // kernel matrix, since it is symmetric. This helps minimize the number of
        // kernel evaluations.
#pragma omp parallel for schedule(dynamic)
        for (size_t i = 0; i < data.n_cols; ++i)
        {
            for (size_t j = i; j < data.n_cols; ++j)
            {
                // Evaluate the kernel on these two points.
                kernelMatrix.at(i, j) = kernel.Evaluate(data.unsafe_col(i), data.unsafe_col(j));
            }
        }

        // Copy to the lower triangular part of the matrix.
        for (size_t i = 1; i < data.n_cols; ++i)
            for (size_t j = 0; j < i; ++j)
                kernelMatrix.at(i, j) = kernelMatrix.at(j, i);

        return kernelMatrix;
    }
};
}  // namespace chenx
