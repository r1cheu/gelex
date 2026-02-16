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

#ifndef GELEX_UTILS_RUNNING_STATS_H_
#define GELEX_UTILS_RUNNING_STATS_H_

#include <cstddef>
#include <type_traits>

#include <Eigen/Core>

#include "gelex/exception.h"

namespace gelex
{

struct RunningStatsResult
{
    Eigen::VectorXd mean;
    Eigen::VectorXd stddev;
};

class RunningStats
{
   public:
    RunningStats() = default;

    template <typename Derived>
        requires std::is_arithmetic_v<typename Derived::Scalar>
    auto update(const Eigen::DenseBase<Derived>& block) -> void
    {
        if (rows_ != 0 && block.rows() != rows_)
        {
            throw InvalidInputException(
                "Row size mismatch in RunningStats::update");
        }

        if (block.cols() == 0)
        {
            return;
        }

        if (rows_ == 0 && block.rows() == 0)
        {
            throw InvalidInputException(
                "Zero-row block in RunningStats::update");
        }

        if (!block.derived().allFinite())
        {
            throw InvalidInputException(
                "Non-finite value in RunningStats::update");
        }

        if (rows_ == 0)
        {
            rows_ = block.rows();
            mean_ = Eigen::VectorXd::Zero(rows_);
            m2_ = Eigen::VectorXd::Zero(rows_);
        }

        if (column_buffer_.size() != rows_)
        {
            column_buffer_.resize(rows_);
            delta_buffer_.resize(rows_);
            delta2_buffer_.resize(rows_);
        }

        for (Eigen::Index col = 0; col < block.cols(); ++col)
        {
            ++count_;
            auto inv_count = 1.0 / static_cast<double>(count_);

            column_buffer_ = block.derived().col(col).template cast<double>();
            delta_buffer_ = column_buffer_ - mean_;
            mean_ += delta_buffer_ * inv_count;
            delta2_buffer_ = column_buffer_ - mean_;
            m2_ += delta_buffer_.cwiseProduct(delta2_buffer_);
        }
    }

    auto result() const -> RunningStatsResult;

   private:
    Eigen::Index rows_{0};
    std::size_t count_{0};
    Eigen::VectorXd mean_;
    Eigen::VectorXd m2_;
    Eigen::VectorXd column_buffer_;
    Eigen::VectorXd delta_buffer_;
    Eigen::VectorXd delta2_buffer_;
};

}  // namespace gelex

#endif  // GELEX_UTILS_RUNNING_STATS_H_
