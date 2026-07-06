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

#ifndef GELEX_INFRA_STATS_NORMAL_SAMPLER_H_
#define GELEX_INFRA_STATS_NORMAL_SAMPLER_H_

#include <cassert>
#include <cmath>
#include <concepts>
#include <random>

namespace gelex
{

template <std::floating_point T>
class NormalSampler
{
   public:
    using Scalar = T;

    struct Params
    {
        T mean{};
        T var{};
    };

    struct Kernel
    {
        T quadratic{};
        T linear{};
        T scale{};

        static auto from_params(const Params& params) -> Kernel
        {
            return {T{1}, params.mean, params.var};
        }
    };

    struct Posterior
    {
        Params params{};
        T log_likelihood_kernel{};
    };

    explicit NormalSampler(T prior_var) : prior_var_(prior_var)
    {
        assert(
            (prior_var_ >= T{0})
            && "NormalSampler: prior variance must be non-negative");
    }

    auto set_prior_var(T value) -> NormalSampler&
    {
        assert(
            (value >= T{0})
            && "NormalSampler: prior variance must be non-negative");
        prior_var_ = value;
        return *this;
    }

    auto posterior(const Params& likelihood) const -> Params
    {
        return posterior(Kernel::from_params(likelihood));
    }

    auto posterior(const Kernel& likelihood) const -> Params
    {
        assert(
            (prior_var_ >= T{0})
            && "NormalSampler: prior variance must be non-negative");
        assert(
            (likelihood.quadratic >= T{0})
            && "NormalSampler: likelihood quadratic term must be non-negative");
        assert(
            (likelihood.scale > T{0})
            && "NormalSampler: likelihood scale must be positive");

        const T denominator
            = likelihood.scale + (likelihood.quadratic * prior_var_);
        assert(
            (denominator > T{0})
            && "NormalSampler: posterior denominator must be positive");

        const T posterior_factor = prior_var_ / denominator;
        return {
            likelihood.linear * posterior_factor,
            likelihood.scale * posterior_factor,
        };
    }

    auto posterior_with_logL(const Kernel& likelihood) const -> Posterior
    {
        const auto posterior_params = posterior(likelihood);
        const T log_likelihood_kernel
            = -T{0.5}
              * (std::log(
                     ((likelihood.quadratic * prior_var_) / likelihood.scale)
                     + T{1})
                 - ((posterior_params.mean * likelihood.linear)
                    / likelihood.scale));
        return {posterior_params, log_likelihood_kernel};
    }

    auto draw(const Params& posterior, std::mt19937_64& rng) -> T
    {
        assert(
            (posterior.var >= T{0})
            && "NormalSampler: posterior variance must be non-negative");
        return (normal_(rng) * std::sqrt(posterior.var)) + posterior.mean;
    }

    auto operator()(const Params& likelihood, std::mt19937_64& rng) -> T
    {
        return draw(posterior(likelihood), rng);
    }

    auto operator()(const Kernel& likelihood, std::mt19937_64& rng) -> T
    {
        return draw(posterior(likelihood), rng);
    }

    auto prior_var() const -> T { return prior_var_; }

    auto reset() -> void { normal_.reset(); }

   private:
    T prior_var_;
    std::normal_distribution<T> normal_{T{0}, T{1}};
};

}  // namespace gelex

#endif  // GELEX_INFRA_STATS_NORMAL_SAMPLER_H_
