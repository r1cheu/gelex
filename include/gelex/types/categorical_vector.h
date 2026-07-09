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

#ifndef GELEX_TYPES_CATEGORICAL_VECTOR_H_
#define GELEX_TYPES_CATEGORICAL_VECTOR_H_

#include <Eigen/Core>
#include <utility>

#include "gelex/exception.h"

namespace gelex
{

class CategoricalVector
{
   public:
    CategoricalVector() = default;
    CategoricalVector(Eigen::Index size, Eigen::Index category_count)
        : values_(Eigen::VectorXi::Zero(size)), category_count_(category_count)
    {
        if (size < 0)
        {
            throw GelexException(
                "CategoricalVector: size must be non-negative");
        }
        if (category_count_ <= 0)
        {
            throw GelexException(
                "CategoricalVector: category count must be positive");
        }
    }
    CategoricalVector(Eigen::VectorXi values, Eigen::Index category_count)
        : values_(std::move(values)), category_count_(category_count)
    {
        if (category_count_ <= 0)
        {
            throw GelexException(
                "CategoricalVector: category count must be positive");
        }
        if ((values_.array() < 0).any()
            || (values_.array() >= category_count_).any())
        {
            throw GelexException("CategoricalVector: category out of range");
        }
    }

    auto operator=(Eigen::VectorXi new_values) -> CategoricalVector&
    {
        if (category_count_ <= 0)
        {
            throw GelexException(
                "CategoricalVector: category count must be positive");
        }
        if (new_values.rows() != values_.rows()
            || new_values.cols() != values_.cols())
        {
            throw GelexException("CategoricalVector: shape changed");
        }
        if ((new_values.array() < 0).any()
            || (new_values.array() >= category_count_).any())
        {
            throw GelexException("CategoricalVector: category out of range");
        }
        values_ = std::move(new_values);
        return *this;
    }

    [[nodiscard]] auto data() const noexcept -> const int*
    {
        return values_.data();
    }

    [[nodiscard]] auto size() const noexcept -> Eigen::Index
    {
        return values_.size();
    }

    [[nodiscard]] auto rows() const noexcept -> Eigen::Index
    {
        return values_.rows();
    }

    [[nodiscard]] auto cols() const noexcept -> Eigen::Index
    {
        return values_.cols();
    }

    auto operator()(Eigen::Index i) noexcept -> int& { return values_(i); }
    [[nodiscard]] auto operator()(Eigen::Index i) const noexcept -> int
    {
        return values_(i);
    }

    [[nodiscard]] auto values() const noexcept -> const Eigen::VectorXi&
    {
        return values_;
    }

    [[nodiscard]] auto category_count() const noexcept -> Eigen::Index
    {
        return category_count_;
    }

   private:
    Eigen::VectorXi values_;
    Eigen::Index category_count_{};
};

}  // namespace gelex

#endif  // GELEX_TYPES_CATEGORICAL_VECTOR_H_
