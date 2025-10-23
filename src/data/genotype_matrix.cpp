#include "gelex/data/genotype_matrix.h"

namespace gelex
{

auto GenotypeMatrix::create(
    Eigen::MatrixXd&& data,
    std::unordered_set<int64_t>&& mono_set,
    Eigen::VectorXd&& mean,
    Eigen::VectorXd&& variance) -> GenotypeMatrix
{
    return GenotypeMatrix(
        std::move(data),
        std::move(mono_set),
        std::move(mean),
        std::move(variance));
}

}  // namespace gelex
