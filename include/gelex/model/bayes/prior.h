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

#ifndef GELEX_MODEL_BAYES_PRIOR_H_
#define GELEX_MODEL_BAYES_PRIOR_H_

#include <memory>
#include <ranges>
#include <vector>

#include "gelex/model/bayes/genetic_prior.h"
#include "gelex/model/bayes/prior_parameters.h"

namespace gelex::bayes
{

class BayesPrior
{
   public:
    BayesPrior(
        VarianceParameter random,
        std::vector<std::unique_ptr<GeneticPrior>> genetics,
        VarianceParameter residual);

    BayesPrior(const BayesPrior&) = delete;
    BayesPrior(BayesPrior&&) noexcept = default;

    auto operator=(const BayesPrior&) -> BayesPrior& = delete;
    auto operator=(BayesPrior&&) noexcept -> BayesPrior& = default;

    ~BayesPrior() = default;

    auto random() const -> const VarianceParameter& { return random_; }
    auto residual() const -> const VarianceParameter& { return residual_; }

    auto genetics() const -> decltype(auto)
    {
        return genetics_
               | std::views::transform(
                   [](const auto& prior) -> const GeneticPrior&
                   { return *prior; });
    }

   private:
    VarianceParameter random_;
    std::vector<std::unique_ptr<GeneticPrior>> genetics_;
    VarianceParameter residual_;
};

}  // namespace gelex::bayes

#endif  // GELEX_MODEL_BAYES_PRIOR_H_
