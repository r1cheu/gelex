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

#ifndef GELEX_ALGO_INFER_POSTERIOR_SUMMARY_H_
#define GELEX_ALGO_INFER_POSTERIOR_SUMMARY_H_

#include <optional>
#include <string>
#include <vector>

#include <Eigen/Core>

#include "gelex/algo/infer/fixed_samples.h"
#include "gelex/infra/stats/running_stats.h"

namespace gelex
{

using stats::RunningStatsResult;

struct PosteriorSummary
{
    explicit PosteriorSummary(Eigen::Index n_params)
        : mean(Eigen::VectorXd::Zero(n_params)),
          stddev(Eigen::VectorXd::Zero(n_params))
    {
    }
    explicit PosteriorSummary(RunningStatsResult result);
    PosteriorSummary() = default;

    auto size() const -> Eigen::Index { return mean.size(); }

    Eigen::VectorXd mean;
    Eigen::VectorXd stddev;
};

struct FixedSummary
{
    explicit FixedSummary(const FixedSamples& sample)
        : names(sample.names), levels(sample.levels), coeffs(sample.n_coeffs())
    {
    }

    void compute(const FixedSamples& sample);

    template <typename F>
    void for_each_term(const F& fn) const
    {
        Eigen::Index coeff_idx = 0;
        for (size_t group_idx = 0; group_idx < levels.size(); ++group_idx)
        {
            const auto& group_levels = levels[group_idx];
            if (group_levels)
            {
                for (const auto& level : *group_levels)
                {
                    fn(names[group_idx] + "_" + level, coeff_idx);
                    ++coeff_idx;
                }
            }
            else
            {
                fn(names[group_idx], coeff_idx);
                ++coeff_idx;
            }
        }
    }

    std::vector<std::string> names;
    std::vector<std::optional<std::vector<std::string>>> levels;
    PosteriorSummary coeffs;
};

struct RandomSummary
{
    explicit RandomSummary(const RandomSamples& sample)
        : name(sample.name),
          levels(sample.levels),
          coeffs(sample.n_coeffs()),
          variance(1)
    {
    }

    void compute(const RandomSamples& sample);

    template <typename F>
    void for_each_term(const F& fn) const
    {
        Eigen::Index coeff_idx = 0;
        if (levels)
        {
            for (const auto& level : *levels)
            {
                fn(name.empty() ? level : name + "_" + level, coeff_idx);
                ++coeff_idx;
            }
        }
        else
        {
            fn(name, coeff_idx);
            ++coeff_idx;
        }
    }

    std::string name;
    std::optional<std::vector<std::string>> levels;
    PosteriorSummary coeffs;
    PosteriorSummary variance;
};

}  // namespace gelex

#endif  // GELEX_ALGO_INFER_POSTERIOR_SUMMARY_H_
