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

#include "gelex/simulate/effect_sampler.h"

#include <Eigen/Core>
#include <cmath>
#include <cstddef>
#include <random>
#include <string>
#include <vector>

#include "gelex/simulate/sim_types.h"

namespace gelex
{

auto NormalGenerator::operator()(std::mt19937_64& rng) const
    -> std::vector<CausalSnp>
{
    std::vector<CausalSnp> causals;
    Eigen::Index start = 0;
    for (const auto& [count, variance] : effect_sizes_)
    {
        const Eigen::Index end = start + count;
        if (variance > 0.0)
        {
            std::normal_distribution<double> dist(0.0, std::sqrt(variance));
            causals.reserve(causals.size() + static_cast<std::size_t>(count));
            for (Eigen::Index i = start; i < end; ++i)
            {
                causals.push_back(
                    CausalSnp{
                        .id = std::string(ids_[i]),
                        .effect = dist(rng),
                    });
            }
        }
        start = end;
    }
    return causals;
}

}  // namespace gelex
