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

#include "gelex/infra/stats/detail/running_stats.h"

#include <algorithm>
#include <cmath>
#include <utility>

namespace gelex::detail
{

auto ScalarRunningStats::result() const noexcept -> ScalarRunningStatsResult
{
    assert(count_ > 0);

    ScalarRunningStatsResult output{.mean = mean_, .stddev = 0.0};
    if (count_ > 1)
    {
        output.stddev
            = std::sqrt(std::max(0.0, m2_ / static_cast<double>(count_ - 1)));
    }
    return output;
}

VectorRunningStats::VectorRunningStats(Eigen::Index size)
{
    if (size <= 0)
    {
        throw GelexException("VectorRunningStats requires a positive size");
    }
    mean_ = Eigen::VectorXd::Zero(size);
    m2_ = Eigen::VectorXd::Zero(size);
    delta_.resize(size);
}

auto VectorRunningStats::result() const -> VectorRunningStatsResult
{
    assert(count_ > 0);

    VectorRunningStatsResult output{
        .mean = mean_, .stddev = Eigen::VectorXd::Zero(mean_.size())};
    if (count_ > 1)
    {
        const Eigen::VectorXd variance
            = (m2_ / static_cast<double>(count_ - 1)).cwiseMax(0.0);
        output.stddev = variance.array().sqrt();
    }
    return output;
}

auto RunningStats::result() const -> RunningStatsResult
{
    RunningStatsResult output;

    if (rows_ == 0)
    {
        return output;
    }

    output.mean = mean_;
    output.stddev = Eigen::VectorXd::Zero(rows_);

    if (count_ <= 1)
    {
        return output;
    }

    Eigen::VectorXd variance = m2_ / static_cast<double>(count_ - 1);
    variance = variance.cwiseMax(0.0);
    output.stddev = variance.array().sqrt();
    return output;
}

CategoricalFrequency::CategoricalFrequency(
    Eigen::Index n_items,
    Eigen::Index n_categories)
    : probabilities_(Eigen::MatrixXd::Zero(n_items, n_categories))
{
}

auto CategoricalFrequency::update(
    const Eigen::Ref<const Eigen::VectorXi>& categories) -> void
{
    const auto old_count = count_;
    ++count_;
    probabilities_
        *= static_cast<double>(old_count) / static_cast<double>(count_);

    const double increment = 1.0 / static_cast<double>(count_);
    for (Eigen::Index i = 0; i < categories.size(); ++i)
    {
        probabilities_(i, categories(i)) += increment;
    }
}

auto CategoricalFrequency::take_probabilities() && -> CategoryProbResult
{
    return CategoryProbResult{std::move(probabilities_)};
}

}  // namespace gelex::detail
