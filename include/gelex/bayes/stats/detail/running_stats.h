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

#ifndef GELEX_BAYES_STATS_DETAIL_RUNNING_STATS_H_
#define GELEX_BAYES_STATS_DETAIL_RUNNING_STATS_H_

#include <Eigen/Core>
#include <cassert>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <ranges>
#include <span>
#include <utility>

#include "gelex/bayes/stats/result.h"
#include "gelex/exception.h"

namespace gelex::detail
{

class ScalarRunningStats
{
   public:
    auto update(double value) noexcept -> void
    {
        assert(std::isfinite(value));
        ++count_;
        const double delta = value - mean_;
        mean_ += delta / static_cast<double>(count_);
        const double delta2 = value - mean_;
        m2_ += delta * delta2;
    }

    [[nodiscard]] auto empty() const noexcept -> bool { return count_ == 0; }

    auto result() const noexcept -> ScalarRunningStatsResult;

   private:
    std::size_t count_{0};
    double mean_{0.0};
    double m2_{0.0};
};

class VectorRunningStats
{
   public:
    explicit VectorRunningStats(Eigen::Index size);

    auto update(const Eigen::Ref<const Eigen::VectorXd>& value) noexcept -> void
    {
        assert(value.size() == mean_.size());
        assert(value.allFinite());

        ++count_;
        const double inv_count = 1.0 / static_cast<double>(count_);
        delta_ = value - mean_;
        mean_ += delta_ * inv_count;
        m2_.array() += delta_.array() * (value.array() - mean_.array());
    }

    [[nodiscard]] auto empty() const noexcept -> bool { return count_ == 0; }

    auto result() const -> VectorRunningStatsResult;

   private:
    std::size_t count_{0};
    Eigen::VectorXd mean_;
    Eigen::VectorXd m2_;
    Eigen::VectorXd delta_;
};

template <std::size_t CategoryCount>
    requires(
        CategoryCount > 1
        && CategoryCount <= static_cast<std::size_t>(
                                std::numeric_limits<std::uint8_t>::max())
                                + 1)
class CategoryRunningStats
{
    using CountMatrix = Eigen::Matrix<
        std::uint32_t,
        Eigen::Dynamic,
        static_cast<int>(CategoryCount),
        Eigen::RowMajor>;

   public:
    explicit CategoryRunningStats(Eigen::Index size)
    {
        if (size <= 0)
        {
            throw GelexException(
                "CategoryRunningStats requires a positive size");
        }
        counts_
            = CountMatrix::Zero(size, static_cast<Eigen::Index>(CategoryCount));
    }

    auto update(std::span<const std::uint8_t> values) noexcept -> void
    {
        assert(std::cmp_equal(values.size(), counts_.rows()));
        assert(count_ < std::numeric_limits<std::uint32_t>::max());

        ++count_;
        for (const auto [index, category] : values | std::views::enumerate)
        {
            assert(static_cast<std::size_t>(category) < CategoryCount);
            ++counts_(
                static_cast<Eigen::Index>(index),
                static_cast<Eigen::Index>(category));
        }
    }

    [[nodiscard]] auto empty() const noexcept -> bool { return count_ == 0; }

    auto result() const -> CategoryRunningStatsResult
    {
        assert(count_ > 0);

        Eigen::MatrixXd probabilities = counts_.template cast<double>();
        probabilities /= static_cast<double>(count_);
        return CategoryRunningStatsResult{
            .probabilities = std::move(probabilities)};
    }

   private:
    CountMatrix counts_;
    std::uint32_t count_{0};
};

}  // namespace gelex::detail

#endif  // GELEX_BAYES_STATS_DETAIL_RUNNING_STATS_H_
