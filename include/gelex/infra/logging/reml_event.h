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

#ifndef GELEX_INFRA_LOGGING_REML_EVENT_H_
#define GELEX_INFRA_LOGGING_REML_EVENT_H_

#include <cstddef>
#include <functional>
#include <string>
#include <variant>
#include <vector>

#include "gelex/types/freq_effect.h"

namespace gelex
{

class FreqModel;
class FreqState;

struct VarianceComponent
{
    freq::GrmType type{freq::GrmType::Unknown};
    double variance{};
    double heritability{};
};

struct LocoRemlResult
{
    std::string chr_name;
    double loglike{};
    std::vector<VarianceComponent> genetic;
    double residual_variance{};
    bool converged{true};

    auto total_h2() const -> double
    {
        double sum = 0.0;
        for (const auto& g : genetic)
        {
            sum += g.heritability;
        }
        return sum;
    }
};

struct RemlEmInitEvent
{
    double loglike;
    std::vector<double> init_variances;
};

struct RemlIterationEvent
{
    size_t iter;
    double loglike;
    std::vector<std::string> labels;
    std::vector<double> variances;
};

struct RemlCompleteEvent
{
    const FreqModel* model;
    const FreqState* state;
    bool converged;
    size_t iter_count;
    size_t max_iter;
    double loglike;
};

using RemlEvent
    = std::variant<RemlEmInitEvent, RemlIterationEvent, RemlCompleteEvent>;

using RemlObserver = std::function<void(const RemlEvent&)>;

}  // namespace gelex

#endif  // GELEX_INFRA_LOGGING_REML_EVENT_H_
