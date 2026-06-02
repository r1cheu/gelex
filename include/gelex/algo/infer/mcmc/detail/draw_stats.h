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

#ifndef GELEX_ALGO_INFER_MCMC_DETAIL_DRAW_STATS_H_
#define GELEX_ALGO_INFER_MCMC_DETAIL_DRAW_STATS_H_

#include <utility>

#include <Eigen/Core>

#include "gelex/infra/stats/detail/running_stats.h"
#include "gelex/infra/stats/result.h"

namespace gelex::mcmc::detail
{

struct ScalarDrawStats
{
    auto store(double value) -> void { stats_.update(value); }
    auto stats() const -> gelex::stats::RunningStatsResult
    {
        return stats_.result();
    }

   private:
    gelex::stats::detail::RunningStats stats_;
};

struct VectorDrawStats
{
    auto store(const Eigen::Ref<const Eigen::VectorXd>& value) -> void
    {
        stats_.update(value);
    }
    auto stats() const -> gelex::stats::RunningStatsResult
    {
        return stats_.result();
    }

   private:
    gelex::stats::detail::RunningStats stats_;
};

struct CategoricalDrawStats
{
    CategoricalDrawStats(Eigen::Index n_items, Eigen::Index n_categories)
        : frequency_(n_items, n_categories)
    {
    }

    auto store(const Eigen::Ref<const Eigen::VectorXi>& categories) -> void
    {
        frequency_.update(categories);
    }
    auto probabilities() const& -> const Eigen::MatrixXd&
    {
        return frequency_.probabilities();
    }
    auto probabilities() && -> Eigen::MatrixXd
    {
        return std::move(frequency_).probabilities();
    }

   private:
    gelex::stats::detail::CategoricalFrequency frequency_;
};

}  // namespace gelex::mcmc::detail

#endif  // GELEX_ALGO_INFER_MCMC_DETAIL_DRAW_STATS_H_
