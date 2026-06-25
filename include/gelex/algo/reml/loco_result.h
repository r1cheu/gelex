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

#ifndef GELEX_ALGO_REML_LOCO_RESULT_H_
#define GELEX_ALGO_REML_LOCO_RESULT_H_

#include <string>
#include <vector>

namespace gelex
{

struct VarianceComponent
{
    std::string name;
    double variance{};
    double variance_se{};
    double variance_ratio{};
    double variance_ratio_se{};
};

struct LocoRemlResult
{
    std::string chr_name;
    double loglike{};
    std::vector<VarianceComponent> random;
    double residual_variance{};
    double residual_variance_se{};
    bool converged{true};

    auto total_ratio() const -> double
    {
        double sum = 0.0;
        for (const auto& r : random)
        {
            sum += r.variance_ratio;
        }
        return sum;
    }
};

}  // namespace gelex

#endif  // GELEX_ALGO_REML_LOCO_RESULT_H_
