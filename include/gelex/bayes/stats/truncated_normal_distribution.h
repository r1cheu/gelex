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

#ifndef GELEX_BAYES_STATS_TRUNCATED_NORMAL_DISTRIBUTION_H_
#define GELEX_BAYES_STATS_TRUNCATED_NORMAL_DISTRIBUTION_H_

#include <cassert>
#include <cmath>
#include <concepts>
#include <cstdint>
#include <limits>
#include <random>

namespace gelex
{

enum class HalfLine : std::uint8_t
{
    Negative,
    Positive,
};

template <std::floating_point T = double>
class TruncatedNormalDistribution
{
   public:
    using result_type = T;

    class param_type
    {
       public:
        using distribution_type = TruncatedNormalDistribution<T>;

        param_type() : param_type(T{0}) {}

        // NOLINTBEGIN(bugprone-easily-swappable-parameters)
        explicit param_type(
            T mean,
            T stddev = T{1},
            HalfLine support = HalfLine::Positive)
            : mean_{mean}, stddev_{stddev}, support_{support}
        {
            assert(stddev_ > T{0});
        }
        // NOLINTEND(bugprone-easily-swappable-parameters)

        [[nodiscard]] auto mean() const -> T { return mean_; }
        [[nodiscard]] auto stddev() const -> T { return stddev_; }
        [[nodiscard]] auto support() const -> HalfLine { return support_; }

        friend auto operator==(const param_type& lhs, const param_type& rhs)
            -> bool = default;

       private:
        friend class TruncatedNormalDistribution<T>;

        T mean_;
        T stddev_;
        HalfLine support_;
    };

    TruncatedNormalDistribution() = default;

    // NOLINTBEGIN(bugprone-easily-swappable-parameters)
    explicit TruncatedNormalDistribution(
        T mean,
        T stddev = T{1},
        HalfLine support = HalfLine::Positive)
        : parameters_{mean, stddev, support}
    {
    }
    // NOLINTEND(bugprone-easily-swappable-parameters)

    explicit TruncatedNormalDistribution(const param_type& parameters)
        : parameters_{parameters}
    {
    }

    auto reset() -> void
    {
        normal_.reset();
        uniform_.reset();
    }

    [[nodiscard]] auto mean() const -> T { return parameters_.mean(); }
    [[nodiscard]] auto stddev() const -> T { return parameters_.stddev(); }
    [[nodiscard]] auto support() const -> HalfLine
    {
        return parameters_.support();
    }

    [[nodiscard]] auto param() const -> param_type { return parameters_; }

    auto param(const param_type& parameters) -> void
    {
        parameters_ = parameters;
    }

    [[nodiscard]] auto min() const -> result_type
    {
        return support() == HalfLine::Positive
                   ? T{0}
                   : std::numeric_limits<T>::lowest();
    }

    [[nodiscard]] auto max() const -> result_type
    {
        return support() == HalfLine::Positive ? std::numeric_limits<T>::max()
                                               : T{0};
    }

    template <typename Generator>
    auto operator()(Generator& rng) -> result_type
    {
        return (*this)(rng, parameters_);
    }

    template <typename Generator>
    auto operator()(Generator& rng, const param_type& parameters) -> result_type
    {
        const T standardized_boundary
            = parameters.support_ == HalfLine::Positive
                  ? -parameters.mean_ / parameters.stddev_
                  : parameters.mean_ / parameters.stddev_;

        if (standardized_boundary <= T{0})
        {
            return draw_rejection(rng, parameters);
        }
        return draw_tail(rng, parameters, standardized_boundary);
    }

   private:
    template <typename Generator>
    auto draw_rejection(Generator& rng, const param_type& parameters) -> T
    {
        while (true)
        {
            const T value
                = parameters.mean_ + (parameters.stddev_ * normal_(rng));
            if (parameters.support_ == HalfLine::Positive && value > T{0})
            {
                return value;
            }
            if (parameters.support_ == HalfLine::Negative && value < T{0})
            {
                return value;
            }
        }
    }

    template <typename Generator>
    auto draw_tail(
        Generator& rng,
        const param_type& parameters,
        T standardized_boundary) -> T
    {
        const T rate
            = (standardized_boundary
               + std::sqrt(
                   (standardized_boundary * standardized_boundary) + T{4}))
              / T{2};
        std::exponential_distribution<T> exponential{rate};

        while (true)
        {
            const T candidate = standardized_boundary + exponential(rng);
            const T difference = candidate - rate;
            const T log_acceptance = -T{0.5} * difference * difference;
            if (std::log(uniform_(rng)) <= log_acceptance)
            {
                const T standardized_value
                    = parameters.support_ == HalfLine::Positive ? candidate
                                                                : -candidate;
                return parameters.mean_
                       + (parameters.stddev_ * standardized_value);
            }
        }
    }

    param_type parameters_;
    std::normal_distribution<T> normal_{T{0}, T{1}};
    std::uniform_real_distribution<T> uniform_{T{0}, T{1}};
};

}  // namespace gelex

#endif  // GELEX_BAYES_STATS_TRUNCATED_NORMAL_DISTRIBUTION_H_
