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

#include "gelex/algo/infer/posterior_summary.h"

#include "gelex/infra/stats/running_stats.h"

namespace gelex
{

PosteriorSummary::PosteriorSummary(RunningStatsResult result)
    : mean(std::move(result.mean)), stddev(std::move(result.stddev))
{
}

void FixedSummary::compute(const FixedSamples& sample)
{
    coeffs = PosteriorSummary(sample.coeffs());
}

void RandomSummary::compute(const RandomSamples& sample)
{
    coeffs = PosteriorSummary(sample.coeffs());
    variance = PosteriorSummary(sample.variance());
}

}  // namespace gelex
