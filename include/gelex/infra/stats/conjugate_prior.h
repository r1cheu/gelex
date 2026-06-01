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

#ifndef GELEX_INFRA_STATS_CONJUGATE_PRIOR_H_
#define GELEX_INFRA_STATS_CONJUGATE_PRIOR_H_

#include <cassert>
#include <cmath>
#include <concepts>
#include <cstdint>
#include <numbers>
#include <random>
#include <type_traits>
#include <utility>

#include <Eigen/Core>

// ---------------------------------------------------------------------------
// Conjugate posterior samplers.
// Templated on floating-point scalar T (float, double, long double).
// Each sampler stores prior hyper-parameters; operator() takes sufficient
// statistics + RNG and returns one draw from the posterior.
// ---------------------------------------------------------------------------

namespace gelex::stats
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
            && "NormalSampler: likelihood quadratic term must be "
               "non-negative");
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

template <std::floating_point T>
class DirichletSampler
{
   public:
    using Scalar = T;
    using Likelihood = Eigen::Ref<const Eigen::VectorXi>;

    explicit DirichletSampler(Eigen::VectorX<T> alpha)
        : alpha_(std::move(alpha))
    {
        assert(
            (alpha_.size() >= 2)
            && "DirichletSampler: alpha must have at least two components");
        assert(
            ((alpha_.array() > T{0}).all())
            && "DirichletSampler: alpha components must be positive");
    }

    auto operator()(const Likelihood& counts, std::mt19937_64& rng)
        -> Eigen::VectorX<T>
    {
        assert(
            (counts.size() == alpha_.size())
            && "DirichletSampler: counts size must match alpha size");
        assert(
            ((counts.array() >= 0).all())
            && "DirichletSampler: counts must be non-negative");

        using ParamT = typename std::gamma_distribution<T>::param_type;
        Eigen::VectorX<T> out(alpha_.size());
        T sum = T{0};
        for (Eigen::Index i = 0; i < alpha_.size(); ++i)
        {
            const T a = alpha_(i) + static_cast<T>(counts(i));
            const T x = gamma_(rng, ParamT{a, T{1}});
            out(i) = x;
            sum += x;
        }
        out /= sum;
        return out;
    }

    auto k() const -> Eigen::Index { return alpha_.size(); }
    auto alpha() const -> const Eigen::VectorX<T>& { return alpha_; }
    auto reset() -> void { gamma_.reset(); }

   private:
    Eigen::VectorX<T> alpha_;
    std::gamma_distribution<T> gamma_{T{1}, T{1}};
};

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

    // Allows improper priors (e.g. nu0 = -2, s2_0 = 0 → reference prior
    // p(σ²) ∝ 1/σ²); posterior validity is checked at draw time via nu1 > 0.
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
        T log_tail_pos;
        T log_tail_neg;
    };

    explicit HalfNormalSampler(T prior_var) : normal_(prior_var) {}

    auto set_prior_var(T value) -> HalfNormalSampler&
    {
        normal_.set_prior_var(value);
        return *this;
    }

    auto posterior(const Kernel& likelihood) const -> Params
    {
        return normal_.posterior(likelihood);
    }

    auto posterior_with_logL(const Kernel& likelihood) const -> Posterior
    {
        const auto [params, log_kernel]
            = normal_.posterior_with_logL(likelihood);
        const T z = params.mean / std::sqrt(params.var);
        return {
            params,
            log_kernel + std::log(T{2}),
            log_phi(z),
            log_phi(-z),
        };
    }

    auto draw(const Params& post, std::int8_t sign, std::mt19937_64& rng) -> T
    {
        assert(
            (sign == 1 || sign == -1)
            && "HalfNormalSampler: sign must be +1 or -1");

        const T mu = post.mean;
        const T sigma = std::sqrt(post.var);
        const T alpha = (sign == 1) ? (-mu / sigma) : (mu / sigma);

        if (alpha <= T{0})
        {
            return draw_rejection(post, sign, rng);
        }
        return draw_devroye(post, alpha, sign, rng);
    }

    auto reset() -> void { normal_.reset(); }

   private:
    NormalSampler<T> normal_;
    std::normal_distribution<T> normal_dist_{T{0}, T{1}};
    std::uniform_real_distribution<T> uniform_{T{0}, T{1}};

    static auto log_phi(T z) -> T
    {
        constexpr T SQRT_2 = std::numbers::sqrt2_v<T>;
        const T erfc_val = std::erfc(-z / SQRT_2);
        if (erfc_val > T{0})
        {
            return std::log(T{0.5}) + std::log(erfc_val);
        }
        return log_phi_asymptotic(z);
    }

    // Mills-ratio expansion to O(1/z^8) for numerical stability when erfc
    // underflows
    static auto log_phi_asymptotic(T z) -> T
    {
        constexpr T LOG_2_PI = std::log(T{2} * std::numbers::pi_v<T>);
        const T z2 = z * z;
        const T z4 = z2 * z2;
        const T z6 = z4 * z2;
        const T correction
            = std::log1p((-T{1} / z2) + (T{3} / z4) - (T{15} / z6));
        return (-T{0.5} * z2) - (T{0.5} * LOG_2_PI) - std::log(-z) + correction;
    }

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

    // Robert (1995) exponential-proposal rejection sampler for truncated N(0,1)
    // on (alpha, inf)
    auto draw_devroye(
        const Params& post,
        T alpha,
        std::int8_t sign,
        std::mt19937_64& rng) -> T
    {
        const T mu = post.mean;
        const T sigma = std::sqrt(post.var);
        // optimal rate: alpha_star = (alpha + sqrt(alpha^2 + 4)) / 2
        const T alpha_star = (alpha + std::sqrt((alpha * alpha) + T{4})) / T{2};
        std::exponential_distribution<T> exp_dist{alpha_star};

        while (true)
        {
            const T e = exp_dist(rng);
            const T z_cand = alpha + e;
            const T diff = z_cand - alpha_star;
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

}  // namespace gelex::stats

#endif  // GELEX_INFRA_STATS_CONJUGATE_PRIOR_H_
