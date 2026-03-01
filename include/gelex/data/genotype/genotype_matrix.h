/*
 * Copyright 2026 RuLei Chen
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 * http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

#ifndef GELEX_DATA_GENOTYPE_MATRIX_H_
#define GELEX_DATA_GENOTYPE_MATRIX_H_

#include <algorithm>
#include <vector>

#include <Eigen/Core>

namespace gelex
{

class GenotypeMatrix
{
   public:
    GenotypeMatrix(
        Eigen::MatrixXd&& data,
        std::vector<int64_t>&& mono_indices,
        Eigen::VectorXd&& mean,
        Eigen::VectorXd&& stddev);

    GenotypeMatrix(const GenotypeMatrix&) = delete;
    GenotypeMatrix(GenotypeMatrix&&) noexcept = default;
    GenotypeMatrix& operator=(const GenotypeMatrix&) = delete;
    GenotypeMatrix& operator=(GenotypeMatrix&&) noexcept = default;
    ~GenotypeMatrix() = default;

    [[nodiscard]] const Eigen::MatrixXd& matrix() const noexcept
    {
        return data_;
    }

    [[nodiscard]] bool is_monomorphic(Eigen::Index marker_idx) const noexcept
    {
        return std::ranges::binary_search(mono_indices_, marker_idx);
    }

    [[nodiscard]] const Eigen::VectorXd& mean() const noexcept { return mean_; }
    [[nodiscard]] const Eigen::VectorXd& stddev() const noexcept
    {
        return stddev_;
    }

    [[nodiscard]] int64_t num_mono() const noexcept
    {
        return static_cast<int64_t>(mono_indices_.size());
    }
    [[nodiscard]] int64_t rows() const noexcept { return data_.rows(); }
    [[nodiscard]] int64_t cols() const noexcept { return data_.cols(); }

   private:
    Eigen::MatrixXd data_;
    std::vector<int64_t> mono_indices_;
    Eigen::VectorXd mean_;
    Eigen::VectorXd stddev_;

    void validate_dimensions() const;
};

}  // namespace gelex

#endif  // GELEX_DATA_GENOTYPE_MATRIX_H_
