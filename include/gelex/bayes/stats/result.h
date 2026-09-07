// Copyright 2026 RuLei Chen
// SPDX-License-Identifier: Apache-2.0

#ifndef GELEX_BAYES_STATS_RESULT_H_
#define GELEX_BAYES_STATS_RESULT_H_

#include <Eigen/Core>

namespace gelex
{

struct ScalarRunningStatsResult
{
    double mean{0.0};
    double stddev{0.0};
};

struct VectorRunningStatsResult
{
    Eigen::VectorXd mean;
    Eigen::VectorXd stddev;
};

struct CategoryRunningStatsResult
{
    Eigen::MatrixXd probabilities;
};

}  // namespace gelex

#endif  // GELEX_BAYES_STATS_RESULT_H_
