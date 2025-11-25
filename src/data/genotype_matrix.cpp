#include "gelex/data/genotype_matrix.h"

namespace gelex
{

GenotypeMatrix::GenotypeMatrix(
    Eigen::MatrixXd&& data,
    std::vector<int64_t>&& mono_indices,
    Eigen::VectorXd&& mean,
    Eigen::VectorXd&& stddev)
    : data_(std::move(data)),
      mono_indices_(std::move(mono_indices)),
      mean_(std::move(mean)),
      stddev_(std::move(stddev))
{
    validate_dimensions();
    std::ranges::sort(mono_indices_);
}

void GenotypeMatrix::validate_dimensions() const
{
    if (data_.cols() != mean_.size())
    {
        throw std::invalid_argument(
            std::format(
                "Dimension mismatch: Matrix cols ({}) != Mean size ({})",
                data_.cols(),
                mean_.size()));
    }

    if (data_.cols() != stddev_.size())
    {
        throw std::invalid_argument(
            std::format(
                "Dimension mismatch: Matrix cols ({}) != Stddev size ({})",
                data_.cols(),
                stddev_.size()));
    }
}

}  // namespace gelex
