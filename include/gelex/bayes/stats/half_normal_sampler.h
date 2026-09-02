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

#ifndef GELEX_BAYES_STATS_HALF_NORMAL_SAMPLER_H_
#define GELEX_BAYES_STATS_HALF_NORMAL_SAMPLER_H_

#include <cassert>
#include <cmath>
#include <concepts>
#include <cstdint>
#include <random>
#include <type_traits>

#include "gelex/bayes/stats/normal_sampler.h"
#include "gelex/infra/normal.h"

namespace gelex
{

template <std::floating_point T>
class HalfNormalSampler
{
   public:
    using Kernel = typename NormalSampler<T>::Kernel;
    using Params = typename NormalSampler<T>::Params;

    static_assert(std::is_same_v<Params, typename NormalSampler<T>::Params>);

    struct Posterior
    {
        Params params;
        T log_marginal_kernel;
        std::int8_t sign;
    };

    explicit HalfNormalSampler(T prior_var) : normal_(prior_var)
    {
        assert(
            (prior_var > T{0})
            && "HalfNormalSampler: prior variance must be positive");
    }

    auto set_prior_var(T value) -> HalfNormalSampler&
    {
        assert(
            (value > T{0})
            && "HalfNormalSampler: prior variance must be positive");
        normal_.set_prior_var(value);
        return *this;
    }

    auto posterior(const Kernel& likelihood) const -> Params
    {
        return normal_.posterior(likelihood);
    }

    auto posterior_with_logL(const Kernel& likelihood, std::int8_t sign) const
        -> Posterior
    {
        assert(
            (sign == 1 || sign == -1)
            && "HalfNormalSampler: sign must be +1 or -1");

        const auto [params, log_kernel]
            = normal_.posterior_with_logL(likelihood);
        assert(
            (params.var > T{0})
            && "HalfNormalSampler: posterior variance must be positive");

        const T z = params.mean / std::sqrt(params.var);
        const T log_tail = (sign == 1) ? log_norm_cdf(z) : log_norm_cdf(-z);
        return {
            params,
            log_kernel + std::log(T{2}) + log_tail,
            sign,
        };
    }

    auto draw(const Params& post, std::int8_t sign, std::mt19937_64& rng) -> T
    {
        assert(
            (sign == 1 || sign == -1)
            && "HalfNormalSampler: sign must be +1 or -1");
        assert(
            (post.var > T{0})
            && "HalfNormalSampler: posterior variance must be positive");

        const T mu = post.mean;
        const T sigma = std::sqrt(post.var);
        const T lower_std = (sign == 1) ? (-mu / sigma) : (mu / sigma);

        if (lower_std <= T{0})
        {
            return draw_rejection(post, sign, rng);
        }
        return draw_devroye(post, lower_std, sign, rng);
    }

    auto draw(const Posterior& post, std::mt19937_64& rng) -> T
    {
        return draw(post.params, post.sign, rng);
    }

    auto reset() -> void
    {
        normal_.reset();
        normal_dist_.reset();
        uniform_.reset();
    }

   private:
    NormalSampler<T> normal_;
    std::normal_distribution<T> normal_dist_{T{0}, T{1}};
    std::uniform_real_distribution<T> uniform_{T{0}, T{1}};

    auto draw_rejection(
        const Params& post,
        std::int8_t sign,
        std::mt19937_64& rng) -> T
    {
        const T mu = post.mean;
        const T sigma = std::sqrt(post.var);
        while (true)
        {
            const T z = normal_dist_(rng);
            const T x = mu + (sigma * z);
            if (sign == 1 && x > T{0})
            {
                return x;
            }
            if (sign == -1 && x < T{0})
            {
                return x;
            }
        }
    }

    auto draw_devroye(
        const Params& post,
        T lower_std,
        std::int8_t sign,
        std::mt19937_64& rng) -> T
    {
        const T mu = post.mean;
        const T sigma = std::sqrt(post.var);
        const T rate
            = (lower_std + std::sqrt((lower_std * lower_std) + T{4})) / T{2};
        std::exponential_distribution<T> exp_dist{rate};

        while (true)
        {
            const T e = exp_dist(rng);
            const T z_cand = lower_std + e;
            const T diff = z_cand - rate;
            const T log_accept = -T{0.5} * (diff * diff);
            const T log_u = std::log(uniform_(rng));
            if (log_u <= log_accept)
            {
                const T z_std = (sign == 1) ? z_cand : -z_cand;
                return mu + (sigma * z_std);
            }
        }
    }
};

}  // namespace gelex

#endif  // GELEX_BAYES_STATS_HALF_NORMAL_SAMPLER_H_
