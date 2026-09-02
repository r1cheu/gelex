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

#ifndef GELEX_BAYES_STATS_DIRICHLET_DISTRIBUTION_H_
#define GELEX_BAYES_STATS_DIRICHLET_DISTRIBUTION_H_

#include <algorithm>
#include <array>
#include <cassert>
#include <cmath>
#include <concepts>
#include <cstddef>
#include <random>
#include <ranges>
#include <utility>

namespace gelex
{

template <std::size_t K, std::floating_point T = double>
    requires(K > 1)
class DirichletDistribution
{
   public:
    using result_type = std::array<T, K>;

    class param_type
    {
       public:
        using distribution_type = DirichletDistribution<K, T>;

        param_type() : concentrations_{} { concentrations_.fill(T{1}); }

        explicit param_type(result_type concentrations)
            : concentrations_{std::move(concentrations)}
        {
            assert(
                std::ranges::all_of(
                    concentrations_,
                    [](T concentration) { return concentration > T{0}; }));
        }

        [[nodiscard]] auto concentrations() const -> result_type
        {
            return concentrations_;
        }

        friend auto operator==(const param_type& lhs, const param_type& rhs)
            -> bool = default;

       private:
        friend class DirichletDistribution<K, T>;

        result_type concentrations_;
    };

    DirichletDistribution() = default;

    explicit DirichletDistribution(result_type concentrations)
        : parameters_{std::move(concentrations)}
    {
    }

    explicit DirichletDistribution(const param_type& parameters)
        : parameters_{parameters}
    {
    }

    auto reset() -> void { gamma_.reset(); }

    [[nodiscard]] auto concentrations() const -> result_type
    {
        return parameters_.concentrations();
    }

    [[nodiscard]] auto param() const -> param_type { return parameters_; }

    auto param(const param_type& parameters) -> void
    {
        parameters_ = parameters;
    }

    template <typename Generator>
    auto operator()(Generator& rng) -> result_type
    {
        return (*this)(rng, parameters_);
    }

    template <typename Generator>
    auto operator()(Generator& rng, const param_type& parameters) -> result_type
    {
        result_type sample{};
        T sum = T{0};
        using GammaParameters = typename std::gamma_distribution<T>::param_type;
        for (auto&& [concentration, value] :
             std::views::zip(parameters.concentrations_, sample))
        {
            value = gamma_(rng, GammaParameters{concentration, T{1}});
            sum += value;
        }
        assert(std::isfinite(sum) && sum > T{0});

        std::ranges::transform(
            sample, sample.begin(), [sum](T value) { return value / sum; });
        return sample;
    }

   private:
    param_type parameters_;
    std::gamma_distribution<T> gamma_{T{1}, T{1}};
};

}  // namespace gelex

#endif  // GELEX_BAYES_STATS_DIRICHLET_DISTRIBUTION_H_
