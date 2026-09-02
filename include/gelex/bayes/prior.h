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

#ifndef GELEX_BAYES_PRIOR_H_
#define GELEX_BAYES_PRIOR_H_

#include <span>
#include <utility>
#include <vector>

#include "gelex/bayes/variance_parameter.h"

namespace gelex
{

template <typename GeneticPrior>
class BayesPrior;

namespace detail
{

template <typename GeneticPrior>
auto make_bayes_prior(
    std::vector<VarianceParameter> random,
    GeneticPrior genetic,
    VarianceParameter residual) -> BayesPrior<GeneticPrior>;

}  // namespace detail

template <typename GeneticPrior>
class BayesPrior
{
   public:
    using genetic_prior_type = GeneticPrior;

    [[nodiscard]] auto random() const noexcept
        -> std::span<const VarianceParameter>
    {
        return random_;
    }

    [[nodiscard]] auto genetic() const noexcept -> const GeneticPrior&
    {
        return genetic_;
    }

    [[nodiscard]] auto residual() const noexcept -> const VarianceParameter&
    {
        return residual_;
    }

   private:
    BayesPrior(
        std::vector<VarianceParameter> random,
        GeneticPrior genetic,
        VarianceParameter residual)
        : random_{std::move(random)},
          genetic_{std::move(genetic)},
          residual_{residual}
    {
    }

    template <typename T>
    friend auto detail::make_bayes_prior(
        std::vector<VarianceParameter> random,
        T genetic,
        VarianceParameter residual) -> BayesPrior<T>;

    std::vector<VarianceParameter> random_;
    GeneticPrior genetic_;
    VarianceParameter residual_;
};

namespace detail
{

template <typename GeneticPrior>
auto make_bayes_prior(
    std::vector<VarianceParameter> random,
    GeneticPrior genetic,
    VarianceParameter residual) -> BayesPrior<GeneticPrior>
{
    return BayesPrior<GeneticPrior>{
        std::move(random), std::move(genetic), residual};
}

}  // namespace detail

}  // namespace gelex

#endif  // GELEX_BAYES_PRIOR_H_
