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

#include <utility>

namespace gelex::stats::detail
{

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

auto CategoricalFrequency::probabilities() const& -> const Eigen::MatrixXd&
{
    return probabilities_;
}

auto CategoricalFrequency::probabilities() && -> Eigen::MatrixXd
{
    return std::move(probabilities_);
}

}  // namespace gelex::stats::detail
