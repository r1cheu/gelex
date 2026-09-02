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

#ifndef GELEX_BAYES_STATS_SCALED_INV_CHI2_DISTRIBUTION_H_
#define GELEX_BAYES_STATS_SCALED_INV_CHI2_DISTRIBUTION_H_

#include <cassert>
#include <concepts>
#include <limits>
#include <random>

namespace gelex
{

template <std::floating_point T = double>
class ScaledInvChi2Distribution
{
   public:
    using result_type = T;

    class param_type
    {
       public:
        using distribution_type = ScaledInvChi2Distribution<T>;

        param_type() : param_type(T{1}, T{1}) {}

        // NOLINTNEXTLINE(bugprone-easily-swappable-parameters)
        param_type(T degrees_of_freedom, T scale)
            : degrees_of_freedom_{degrees_of_freedom}, scale_{scale}
        {
            assert(degrees_of_freedom_ > T{0});
            assert(scale_ > T{0});
        }

        [[nodiscard]] auto degrees_of_freedom() const -> T
        {
            return degrees_of_freedom_;
        }

        [[nodiscard]] auto scale() const -> T { return scale_; }

        friend auto operator==(const param_type& lhs, const param_type& rhs)
            -> bool = default;

       private:
        friend class ScaledInvChi2Distribution<T>;

        T degrees_of_freedom_;
        T scale_;
    };

    ScaledInvChi2Distribution() = default;

    // NOLINTNEXTLINE(bugprone-easily-swappable-parameters)
    ScaledInvChi2Distribution(T degrees_of_freedom, T scale)
        : parameters_{degrees_of_freedom, scale}
    {
    }

    explicit ScaledInvChi2Distribution(const param_type& parameters)
        : parameters_{parameters}
    {
    }

    auto reset() -> void { chi_squared_.reset(); }

    [[nodiscard]] auto degrees_of_freedom() const -> T
    {
        return parameters_.degrees_of_freedom();
    }

    [[nodiscard]] auto scale() const -> T { return parameters_.scale(); }

    [[nodiscard]] auto param() const -> param_type { return parameters_; }

    auto param(const param_type& parameters) -> void
    {
        parameters_ = parameters;
    }

    [[nodiscard]] constexpr auto min() const noexcept -> result_type
    {
        return T{0};
    }

    [[nodiscard]] constexpr auto max() const noexcept -> result_type
    {
        return std::numeric_limits<T>::infinity();
    }

    template <typename Generator>
    auto operator()(Generator& rng) -> result_type
    {
        return (*this)(rng, parameters_);
    }

    template <typename Generator>
    auto operator()(Generator& rng, const param_type& parameters) -> result_type
    {
        using ChiSquaredParameters =
            typename std::chi_squared_distribution<T>::param_type;
        const T degrees_of_freedom = parameters.degrees_of_freedom_;
        return (degrees_of_freedom * parameters.scale_)
               / chi_squared_(rng, ChiSquaredParameters{degrees_of_freedom});
    }

   private:
    param_type parameters_;
    std::chi_squared_distribution<T> chi_squared_{T{1}};
};

}  // namespace gelex

#endif  // GELEX_BAYES_STATS_SCALED_INV_CHI2_DISTRIBUTION_H_
