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

#ifndef GELEX_ALGO_INFER_FIXED_SAMPLES_H_
#define GELEX_ALGO_INFER_FIXED_SAMPLES_H_

#include <optional>
#include <string>
#include <vector>

#include <Eigen/Core>

#include "gelex/infra/stats/running_stats.h"

namespace gelex
{
struct FixedEffect;
}

namespace gelex::bayes
{
struct RandomEffect;
struct FixedState;
struct RandomState;
}  // namespace gelex::bayes

namespace gelex
{

using stats::RunningStats;
using stats::RunningStatsResult;

struct FixedSamples
{
    explicit FixedSamples(const FixedEffect& effect);
    void store(const bayes::FixedState& state);

    auto n_coeffs() const -> Eigen::Index { return n_coeffs_; }
    auto coeffs() const -> RunningStatsResult { return coeffs_stats_.result(); }

    std::vector<std::string> names;
    std::vector<std::optional<std::vector<std::string>>> levels;

   private:
    Eigen::Index n_coeffs_;
    RunningStats coeffs_stats_;
};

struct RandomSamples
{
    explicit RandomSamples(const bayes::RandomEffect& effect);
    void store(const bayes::RandomState& state);

    auto n_coeffs() const -> Eigen::Index { return n_coeffs_; }
    auto coeffs() const -> RunningStatsResult { return coeffs_stats_.result(); }
    auto variance() const -> RunningStatsResult
    {
        return variance_stats_.result();
    }

    std::string name;
    std::optional<std::vector<std::string>> levels;

   private:
    Eigen::Index n_coeffs_;
    RunningStats coeffs_stats_;
    RunningStats variance_stats_;
};

}  // namespace gelex

#endif  // GELEX_ALGO_INFER_FIXED_SAMPLES_H_
