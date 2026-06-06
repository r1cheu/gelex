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

#ifndef GELEX_ALGO_MCMC_DETAIL_RECORDS_H_
#define GELEX_ALGO_MCMC_DETAIL_RECORDS_H_

#include <cstddef>
#include <optional>
#include <string_view>
#include <utility>

#include <Eigen/Core>

#include "gelex/infra/stats/detail/running_stats.h"
#include "gelex/infra/stats/result.h"
#include "gelex/types/categorical_vector.h"

namespace gelex::mcmc
{

class Records;

}  // namespace gelex::mcmc

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

    auto store(const CategoricalVector& categories) -> void
    {
        frequency_.update(categories.values());
    }
    auto take_probabilities() && -> gelex::stats::CategoryProbResult
    {
        return std::move(frequency_).take_probabilities();
    }

   private:
    gelex::stats::detail::CategoricalFrequency frequency_;
};

class ScalarRecord
{
   public:
    ScalarRecord(Records& owner, std::string_view path);

    auto store(Records& owner, double value) -> void;
    auto result() const -> stats::RunningStatsResult;

   private:
    ScalarDrawStats draws_;
    std::optional<std::size_t> draw_handle_;
};

class VectorRecord
{
   public:
    VectorRecord(
        Records& owner,
        std::string_view path,
        const Eigen::Ref<const Eigen::VectorXd>& value);

    auto store(Records& owner, const Eigen::Ref<const Eigen::VectorXd>& value)
        -> void;
    auto result() const -> stats::RunningStatsResult;

   private:
    VectorDrawStats draws_;
    std::optional<std::size_t> draw_handle_;
};

class CategoricalRecord
{
   public:
    CategoricalRecord(
        Records& owner,
        std::string_view path,
        const CategoricalVector& value);

    auto store(Records& owner, const CategoricalVector& value) -> void;
    auto result() && -> stats::CategoryProbResult;

   private:
    CategoricalDrawStats draws_;
    std::optional<std::size_t> draw_handle_;
};

}  // namespace gelex::mcmc::detail

#endif  // GELEX_ALGO_MCMC_DETAIL_RECORDS_H_
