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

#include <Eigen/Core>
#include <cassert>
#include <cstddef>

#include "gelex/infra/stats/result.h"

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

    auto result() const -> VectorRunningStatsResult;

   private:
    std::size_t count_{0};
    Eigen::VectorXd mean_;
    Eigen::VectorXd m2_;
    Eigen::VectorXd delta_;
};

}  // namespace gelex::detail

#endif  // GELEX_INFRA_STATS_DETAIL_RUNNING_STATS_H_
