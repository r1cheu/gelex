// Copyright 2026 RuLei Chen
// SPDX-License-Identifier: Apache-2.0

#include "gelex/bayes/stats/detail/running_stats.h"

#include <algorithm>
#include <cmath>

#include "gelex/exception.h"

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

auto VectorRunningStats::mean_square() const -> Eigen::VectorXd
{
    assert(count_ > 0);
    return mean_.array().square() + (m2_.array() / static_cast<double>(count_));
}

}  // namespace gelex::detail
