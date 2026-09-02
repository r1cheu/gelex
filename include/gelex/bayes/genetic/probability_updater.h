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

#ifndef GELEX_BAYES_GENETIC_PROBABILITY_UPDATER_H_
#define GELEX_BAYES_GENETIC_PROBABILITY_UPDATER_H_

#include <Eigen/Core>
#include <algorithm>
#include <array>
#include <cstddef>
#include <random>

#include "gelex/bayes/genetic/prior.h"
#include "gelex/bayes/stats/beta_sampler.h"
#include "gelex/bayes/stats/dirichlet_sampler.h"

namespace gelex::detail
{

template <UpdatePolicy Policy>
class ProbabilityUpdater;

template <>
class ProbabilityUpdater<UpdatePolicy::Fixed>
{
   public:
    explicit ProbabilityUpdater(
        const ProbabilityParameter<UpdatePolicy::Fixed>& /*parameter*/) noexcept
    {
    }
};

template <>
class ProbabilityUpdater<UpdatePolicy::Sampled>
{
   public:
    explicit ProbabilityUpdater(
        const ProbabilityParameter<UpdatePolicy::Sampled>& parameter)
        : sampler_{parameter.hyperprior.alpha, parameter.hyperprior.beta}
    {
    }

    auto reset() -> void
    {
        sampler_.reset();
        n_success_ = 0;
        n_fail_ = 0;
    }

    auto observe(bool success) noexcept -> void
    {
        if (success)
        {
            ++n_success_;
        }
        else
        {
            ++n_fail_;
        }
    }

    auto update(double& probability, std::mt19937_64& rng) -> void
    {
        probability
            = sampler_({.n_success = n_success_, .n_fail = n_fail_}, rng);
    }

   private:
    BetaSampler<double> sampler_;
    Eigen::Index n_success_{};
    Eigen::Index n_fail_{};
};

template <UpdatePolicy Policy, std::size_t Classes>
class SimplexUpdater;

template <std::size_t Classes>
class SimplexUpdater<UpdatePolicy::Fixed, Classes>
{
   public:
    explicit SimplexUpdater(
        const SimplexParameter<Classes, UpdatePolicy::Fixed>&
        /*parameter*/) noexcept
    {
    }
};

template <std::size_t Classes>
class SimplexUpdater<UpdatePolicy::Sampled, Classes>
{
    static constexpr int class_count = static_cast<int>(Classes);
    using Concentrations = Eigen::Vector<double, class_count>;
    using Counts = Eigen::Vector<int, class_count>;

   public:
    explicit SimplexUpdater(
        const SimplexParameter<Classes, UpdatePolicy::Sampled>& parameter)
        : sampler_{Eigen::Map<const Concentrations>{
              parameter.hyperprior.concentration.data()}}
    {
    }

    auto reset() -> void
    {
        sampler_.reset();
        counts_.setZero();
    }

    auto observe(std::size_t class_index) noexcept -> void
    {
        ++counts_(static_cast<Eigen::Index>(class_index));
    }

    auto update(
        std::array<double, Classes>& probabilities,
        std::mt19937_64& rng) -> void
    {
        const auto sample = sampler_(counts_, rng);
        std::copy_n(sample.data(), Classes, probabilities.begin());
    }

   private:
    DirichletSampler<double> sampler_;
    Counts counts_{Counts::Zero()};
};

}  // namespace gelex::detail

#endif  // GELEX_BAYES_GENETIC_PROBABILITY_UPDATER_H_
