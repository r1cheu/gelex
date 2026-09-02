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

#ifndef GELEX_BAYES_STATS_LOG_CATEGORICAL_DISTRIBUTION_H_
#define GELEX_BAYES_STATS_LOG_CATEGORICAL_DISTRIBUTION_H_

#include <algorithm>
#include <array>
#include <cassert>
#include <cmath>
#include <cstddef>
#include <random>
#include <ranges>
#include <utility>

namespace gelex
{

/**
 * @brief Fixed-size categorical distribution parameterized by unnormalized
 * log weights.
 *
 * `param_type` uses the max-shift from log-sum-exp normalization:
 * @f$ p_k = \exp(L_k - M) / \sum_i \exp(L_i - M) @f$, where
 * @f$ M = \max_i L_i @f$. At least one log weight must be finite;
 * negative infinity represents zero weight.
 */
template <std::size_t K>
    requires(K > 0)
class LogCategoricalDistribution
{
   public:
    using result_type = std::size_t;
    using log_weights_type = std::array<double, K>;
    using probabilities_type = std::array<double, K>;

    class param_type
    {
       public:
        using distribution_type = LogCategoricalDistribution<K>;

        param_type() : param_type(log_weights_type{}) {}

        explicit param_type(log_weights_type log_weights)
            : probabilities_{}, cumulative_probabilities_{}
        {
            const double max_log_weight = std::ranges::max(log_weights);
            assert(std::isfinite(max_log_weight));

            // shift log weights by max_log_weight and exp to get probabilities
            double sum_exp = 0.0;
            for (auto&& [log_weight, probability] :
                 std::views::zip(log_weights, probabilities_))
            {
                probability = std::exp(log_weight - max_log_weight);
                sum_exp += probability;
            }
            assert(std::isfinite(sum_exp) && sum_exp > 0.0);

            double cumulative_probability = 0.0;
            for (auto&& [probability, cumulative] :
                 std::views::zip(probabilities_, cumulative_probabilities_))
            {
                probability /= sum_exp;
                cumulative_probability += probability;
                cumulative = cumulative_probability;
            }
            cumulative_probabilities_.back() = 1.0;
        }

        [[nodiscard]] auto probabilities() const -> probabilities_type
        {
            return probabilities_;
        }

        friend auto operator==(const param_type& lhs, const param_type& rhs)
            -> bool = default;

       private:
        friend class LogCategoricalDistribution<K>;

        probabilities_type probabilities_;
        probabilities_type cumulative_probabilities_;
    };

    LogCategoricalDistribution() = default;

    explicit LogCategoricalDistribution(log_weights_type log_weights)
        : parameters_{std::move(log_weights)}
    {
    }

    explicit LogCategoricalDistribution(const param_type& parameters)
        : parameters_{parameters}
    {
    }

    auto reset() noexcept -> void {}

    [[nodiscard]] auto probabilities() const -> probabilities_type
    {
        return parameters_.probabilities();
    }

    [[nodiscard]] auto param() const -> param_type { return parameters_; }

    auto param(const param_type& parameters) -> void
    {
        parameters_ = parameters;
    }

    [[nodiscard]] constexpr auto min() const noexcept -> result_type
    {
        return 0;
    }

    [[nodiscard]] constexpr auto max() const noexcept -> result_type
    {
        return K - 1;
    }

    template <typename Generator>
    auto operator()(Generator& rng) -> result_type
    {
        return (*this)(rng, parameters_);
    }

    template <typename Generator>
    auto operator()(Generator& rng, const param_type& parameters) -> result_type
    {
        const double threshold
            = std::uniform_real_distribution<double>{0.0, 1.0}(rng);
        for (const auto [class_index, cumulative_probability] :
             parameters.cumulative_probabilities_ | std::views::enumerate)
        {
            if (threshold < cumulative_probability)
            {
                return static_cast<result_type>(class_index);
            }
        }
        return max();
    }

    friend auto operator==(
        const LogCategoricalDistribution& lhs,
        const LogCategoricalDistribution& rhs) -> bool = default;

   private:
    param_type parameters_;
};

template <std::size_t K>
    requires(K > 0)
[[nodiscard]] auto make_log_weights(std::array<double, K> probabilities)
    -> std::array<double, K>
{
    std::ranges::transform(
        probabilities,
        probabilities.begin(),
        [](double probability) { return std::log(probability); });
    return probabilities;
}

template <std::size_t K>
    requires(K > 0)
[[nodiscard]] auto make_mixture_posterior_weights(
    const std::array<double, K>& log_mixture_weights,
    std::array<double, K> component_log_integrals) ->
    typename LogCategoricalDistribution<K>::param_type
{
    for (auto&& [log_integral, log_mixture_weight] :
         std::views::zip(component_log_integrals, log_mixture_weights))
    {
        log_integral += log_mixture_weight;
    }
    return typename LogCategoricalDistribution<K>::param_type{
        std::move(component_log_integrals)};
}

}  // namespace gelex

#endif  // GELEX_BAYES_STATS_LOG_CATEGORICAL_DISTRIBUTION_H_
