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

#ifndef GELEX_INFRA_STATS_BETA_SAMPLER_H_
#define GELEX_INFRA_STATS_BETA_SAMPLER_H_

#include <Eigen/Core>
#include <cassert>
#include <concepts>
#include <random>

namespace gelex
{

template <std::floating_point T>
class BetaSampler
{
   public:
    using Scalar = T;

    struct Likelihood
    {
        Eigen::Index n_success{};
        Eigen::Index n_fail{};
    };

    // NOLINTNEXTLINE(bugprone-easily-swappable-parameters)
    BetaSampler(T alpha, T beta) : alpha_(alpha), beta_(beta)
    {
        assert((alpha_ > T{0}) && "BetaSampler: alpha must be positive");
        assert((beta_ > T{0}) && "BetaSampler: beta must be positive");
    }

    auto operator()(const Likelihood& likelihood, std::mt19937_64& rng) -> T
    {
        assert(
            (likelihood.n_success >= 0)
            && "BetaSampler: success count must be non-negative");
        assert(
            (likelihood.n_fail >= 0)
            && "BetaSampler: failure count must be non-negative");

        const T a = alpha_ + static_cast<T>(likelihood.n_success);
        const T b = beta_ + static_cast<T>(likelihood.n_fail);
        using ParamT = typename std::gamma_distribution<T>::param_type;
        const T x = gamma_(rng, ParamT{a, T{1}});
        const T y = gamma_(rng, ParamT{b, T{1}});
        return x / (x + y);
    }

    auto alpha() const -> T { return alpha_; }
    auto beta() const -> T { return beta_; }
    auto reset() -> void { gamma_.reset(); }

   private:
    T alpha_;
    T beta_;
    std::gamma_distribution<T> gamma_{T{1}, T{1}};
};

}  // namespace gelex

#endif  // GELEX_INFRA_STATS_BETA_SAMPLER_H_
