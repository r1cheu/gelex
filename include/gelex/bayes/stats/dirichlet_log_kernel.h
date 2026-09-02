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

#ifndef GELEX_BAYES_STATS_DIRICHLET_LOG_KERNEL_H_
#define GELEX_BAYES_STATS_DIRICHLET_LOG_KERNEL_H_

#include <algorithm>
#include <array>
#include <cassert>
#include <cstddef>
#include <ranges>
#include <utility>

#include "gelex/bayes/stats/dirichlet_distribution.h"

namespace gelex
{

template <std::size_t K>
    requires(K > 1)
class DirichletLogKernel
{
   public:
    explicit DirichletLogKernel(std::array<double, K> exponents)
        : exponents_{std::move(exponents)}
    {
        assert(
            std::ranges::all_of(
                exponents_, [](double exponent) { return exponent > -1.0; }));
    }

    [[nodiscard]] auto dirichlet_parameters() const ->
        typename DirichletDistribution<K>::param_type
    {
        auto parameters = exponents_;
        // to concentrations, we need to add 1 to each exponent
        std::ranges::transform(
            parameters,
            parameters.begin(),
            [](double exponent) { return exponent + 1.0; });
        return typename DirichletDistribution<K>::param_type{
            std::move(parameters)};
    }

    [[nodiscard]] auto operator+(const DirichletLogKernel& other) const
        -> DirichletLogKernel
    {
        auto exponents = exponents_;
        for (auto&& [exponent, other_exponent] :
             std::views::zip(exponents, other.exponents_))
        {
            exponent += other_exponent;
        }
        return DirichletLogKernel{std::move(exponents)};
    }

   private:
    std::array<double, K> exponents_;
};

template <std::size_t K>
    requires(K > 1)
[[nodiscard]] auto make_categorical_likelihood(
    const std::array<std::size_t, K>& counts) -> DirichletLogKernel<K>
{
    std::array<double, K> exponents{};
    std::ranges::transform(
        counts,
        exponents.begin(),
        [](std::size_t count) { return static_cast<double>(count); });
    return DirichletLogKernel<K>{std::move(exponents)};
}

template <std::size_t K>
    requires(K > 1)
[[nodiscard]] auto make_dirichlet_prior(std::array<double, K> concentrations)
    -> DirichletLogKernel<K>
{
    std::ranges::transform(
        concentrations,
        concentrations.begin(),
        [](double concentration) { return concentration - 1.0; });
    return DirichletLogKernel<K>{std::move(concentrations)};
}

template <std::size_t K>
    requires(K > 1)
[[nodiscard]] auto make_uniform_dirichlet_prior() -> DirichletLogKernel<K>
{
    std::array<double, K> concentrations{};
    concentrations.fill(1.0);
    return make_dirichlet_prior(std::move(concentrations));
}

[[nodiscard]] inline auto make_beta_prior(double alpha, double beta)
    -> DirichletLogKernel<2>
{
    return DirichletLogKernel<2>{{beta - 1.0, alpha - 1.0}};
}

}  // namespace gelex

#endif  // GELEX_BAYES_STATS_DIRICHLET_LOG_KERNEL_H_
