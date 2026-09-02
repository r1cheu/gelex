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

#ifndef GELEX_BAYES_STATS_SCALED_INV_CHI2_SAMPLER_H_
#define GELEX_BAYES_STATS_SCALED_INV_CHI2_SAMPLER_H_

#include <Eigen/Core>
#include <cassert>
#include <concepts>
#include <random>

namespace gelex
{

template <std::floating_point T>
class ScaledInvChi2Sampler
{
   public:
    using Scalar = T;

    struct Likelihood
    {
        Eigen::Index n{};
        T sum_squares{};
    };

    // NOLINTNEXTLINE(bugprone-easily-swappable-parameters)
    ScaledInvChi2Sampler(T nu0, T s2_0) : nu0_(nu0), s2_0_(s2_0)
    {
        assert(
            (s2_0_ >= T{0})
            && "ScaledInvChi2Sampler: prior scale must be non-negative");
    }

    template <typename Prior>
        requires requires(const Prior& prior) {
            { prior.degrees_of_freedom() } -> std::convertible_to<T>;
            { prior.scale() } -> std::convertible_to<T>;
        }
    explicit ScaledInvChi2Sampler(const Prior& prior)
        : ScaledInvChi2Sampler(
              static_cast<T>(prior.degrees_of_freedom()),
              static_cast<T>(prior.scale()))
    {
    }

    auto operator()(const Likelihood& likelihood, std::mt19937_64& rng) -> T
    {
        assert(
            (likelihood.n >= 0)
            && "ScaledInvChi2Sampler: observation count must be non-negative");
        assert(
            (likelihood.sum_squares >= T{0})
            && "ScaledInvChi2Sampler: sum of squares must be non-negative");

        const T nu1 = nu0_ + static_cast<T>(likelihood.n);
        assert(
            (nu1 > T{0})
            && "ScaledInvChi2Sampler: posterior nu must be positive");

        const T s2_1 = ((nu0_ * s2_0_) + likelihood.sum_squares) / nu1;
        assert(
            (s2_1 > T{0})
            && "ScaledInvChi2Sampler: posterior scale must be positive");

        using ParamT = typename std::chi_squared_distribution<T>::param_type;
        return (nu1 * s2_1) / chisq_(rng, ParamT{nu1});
    }

    auto nu0() const -> T { return nu0_; }
    auto s2_0() const -> T { return s2_0_; }
    auto reset() -> void { chisq_.reset(); }

   private:
    T nu0_;
    T s2_0_;
    std::chi_squared_distribution<T> chisq_{T{1}};
};

}  // namespace gelex

#endif  // GELEX_BAYES_STATS_SCALED_INV_CHI2_SAMPLER_H_
