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

#include <concepts>
#include <random>
#include <utility>

#include <Eigen/Core>

#include "gelex/exception.h"

namespace gelex
{
namespace detail
{

struct ScaledInvChiSqParams
{
    double nu{};
    double s2{};
};

struct NormalParams
{
    double mean{};
    double var{};
};

class ScaledInvChiSq
{
   public:
    explicit ScaledInvChiSq(const ScaledInvChiSqParams& prior_params);
    ScaledInvChiSq(double initial_nu, double initial_s2);

    void compute(double sum_of_squared_errors, Eigen::Index num_observations);

    void compute(double single_observation_squared_error);

    double operator()(std::mt19937_64& rng) const;
    auto expected_value() const -> double;
    auto posterior_stddev() const -> double;
    const ScaledInvChiSqParams& prior() { return prior_; }
    const ScaledInvChiSqParams& posterior() { return posterior_; }

   private:
    ScaledInvChiSqParams prior_;
    ScaledInvChiSqParams posterior_;
};

}  // namespace detail

// ---------------------------------------------------------------------------
// Stateless conjugate posterior samplers.
// Templated on floating-point scalar T (float, double, long double).
// Each sampler stores prior hyperparameters; operator() takes sufficient
// statistics + RNG and returns one draw from the posterior. Samplers are
// thread-safe (const) — share one instance across threads, give each thread
// its own RNG.
// ---------------------------------------------------------------------------

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
        if (!(alpha_ > T{0}) || !(beta_ > T{0}))
        {
            throw GelexException(
                "BetaSampler: alpha and beta must be positive");
        }
    }

    auto operator()(const Likelihood& lik, std::mt19937_64& rng) const -> T
    {
        if (lik.n_success < 0 || lik.n_fail < 0)
        {
            throw GelexException("BetaSampler: counts must be non-negative");
        }
        const T a = alpha_ + static_cast<T>(lik.n_success);
        const T b = beta_ + static_cast<T>(lik.n_fail);
        std::gamma_distribution<T> ga{a, T{1}};
        std::gamma_distribution<T> gb{b, T{1}};
        const T x = ga(rng);
        const T y = gb(rng);
        return x / (x + y);
    }

    auto alpha() const -> T { return alpha_; }
    auto beta() const -> T { return beta_; }

   private:
    T alpha_;
    T beta_;
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
        if (alpha_.size() < 2)
        {
            throw GelexException("DirichletSampler: alpha must have size >= 2");
        }
        if ((alpha_.array() <= T{0}).any())
        {
            throw GelexException(
                "DirichletSampler: all alpha components must be positive");
        }
    }

    auto operator()(const Likelihood& counts, std::mt19937_64& rng) const
        -> Eigen::VectorX<T>
    {
        if (counts.size() != alpha_.size())
        {
            throw GelexException(
                "DirichletSampler: counts size must match alpha size");
        }
        if ((counts.array() < 0).any())
        {
            throw GelexException(
                "DirichletSampler: counts must be non-negative");
        }

        Eigen::VectorX<T> out(alpha_.size());
        T sum = T{0};
        for (Eigen::Index i = 0; i < alpha_.size(); ++i)
        {
            const T a = alpha_(i) + static_cast<T>(counts(i));
            std::gamma_distribution<T> g{a, T{1}};
            const T x = g(rng);
            out(i) = x;
            sum += x;
        }
        out /= sum;
        return out;
    }

    auto k() const -> Eigen::Index { return alpha_.size(); }
    auto alpha() const -> const Eigen::VectorX<T>& { return alpha_; }

   private:
    Eigen::VectorX<T> alpha_;
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
        if (s2_0_ < T{0})
        {
            throw GelexException(
                "ScaledInvChi2Sampler: s2_0 must be non-negative");
        }
    }

    auto operator()(const Likelihood& lik, std::mt19937_64& rng) const -> T
    {
        if (lik.n < 0)
        {
            throw GelexException(
                "ScaledInvChi2Sampler: n must be non-negative");
        }
        if (lik.sum_squares < T{0})
        {
            throw GelexException(
                "ScaledInvChi2Sampler: sum_squares must be non-negative");
        }
        const T nu1 = nu0_ + static_cast<T>(lik.n);
        if (!(nu1 > T{0}))
        {
            throw GelexException(
                "ScaledInvChi2Sampler: posterior nu must be positive (improper "
                "prior with n == 0)");
        }
        const T s2_1 = ((nu0_ * s2_0_) + lik.sum_squares) / nu1;
        if (!(s2_1 > T{0}))
        {
            throw GelexException(
                "ScaledInvChi2Sampler: posterior s2 must be positive");
        }
        std::chi_squared_distribution<T> chisq{nu1};
        return (nu1 * s2_1) / chisq(rng);
    }

    auto nu0() const -> T { return nu0_; }
    auto s2_0() const -> T { return s2_0_; }

   private:
    T nu0_;
    T s2_0_;
};

}  // namespace gelex

#endif  // GELEX_INFRA_STATS_CONJUGATE_PRIOR_H_
