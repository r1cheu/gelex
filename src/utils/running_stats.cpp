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

#include "gelex/utils/running_stats.h"

#include <cmath>

namespace gelex
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

}  // namespace gelex
