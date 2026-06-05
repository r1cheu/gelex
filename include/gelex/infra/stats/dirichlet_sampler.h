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

#ifndef GELEX_INFRA_STATS_DIRICHLET_SAMPLER_H_
#define GELEX_INFRA_STATS_DIRICHLET_SAMPLER_H_

#include <cassert>
#include <concepts>
#include <random>
#include <utility>

#include <Eigen/Core>

namespace gelex::stats
{

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

}  // namespace gelex::stats

#endif  // GELEX_INFRA_STATS_DIRICHLET_SAMPLER_H_
