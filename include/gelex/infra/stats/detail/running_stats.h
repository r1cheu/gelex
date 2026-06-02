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

#ifndef GELEX_INFRA_STATS_DETAIL_RUNNING_STATS_H_
#define GELEX_INFRA_STATS_DETAIL_RUNNING_STATS_H_

#include <cmath>
#include <cstddef>
#include <type_traits>

#include <Eigen/Core>

#include "gelex/exception.h"
#include "gelex/infra/stats/result.h"

namespace gelex::stats::detail
{

class RunningStats
{
   public:
    RunningStats() = default;

    auto update(double value) -> void
    {
        if (rows_ == 0)
        {
            rows_ = 1;
            mean_ = Eigen::VectorXd::Zero(1);
            m2_ = Eigen::VectorXd::Zero(1);
        }
        if (rows_ != 1)
        {
            throw GelexException("Scalar update on multi-row RunningStats");
        }
        if (!std::isfinite(value))
        {
            throw GelexException("Non-finite value in RunningStats::update");
        }
        ++count_;
        const double delta = value - mean_(0);
        mean_(0) += delta / static_cast<double>(count_);
        const double delta2 = value - mean_(0);
        m2_(0) += delta * delta2;
    }

    template <typename Derived>
        requires std::is_arithmetic_v<typename Derived::Scalar>
    auto update(const Eigen::DenseBase<Derived>& block) -> void
    {
        if (rows_ != 0 && block.rows() != rows_)
        {
            throw GelexException("Row size mismatch in RunningStats::update");
        }

        if (block.cols() == 0)
        {
            return;
        }

        if (rows_ == 0 && block.rows() == 0)
        {
            throw GelexException("Zero-row block in RunningStats::update");
        }

        if (!block.derived().allFinite())
        {
            throw GelexException("Non-finite value in RunningStats::update");
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

class CategoricalFrequency
{
   public:
    CategoricalFrequency(Eigen::Index n_items, Eigen::Index n_categories);

    auto update(const Eigen::Ref<const Eigen::VectorXi>& categories) -> void;
    auto probabilities() const& -> const Eigen::MatrixXd&;
    auto probabilities() && -> Eigen::MatrixXd;

   private:
    Eigen::MatrixXd probabilities_;
    std::size_t count_{0};
};

}  // namespace gelex::stats::detail

#endif  // GELEX_INFRA_STATS_DETAIL_RUNNING_STATS_H_
