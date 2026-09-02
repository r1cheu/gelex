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

#ifndef GELEX_BAYES_GENETIC_KERNEL_ALLOCATION_UPDATER_H_
#define GELEX_BAYES_GENETIC_KERNEL_ALLOCATION_UPDATER_H_

#include <Eigen/Core>
#include <algorithm>
#include <array>
#include <cassert>
#include <cmath>
#include <cstddef>
#include <random>

#include "gelex/bayes/genetic/parameter.h"
#include "gelex/bayes/stats/beta_sampler.h"
#include "gelex/bayes/stats/dirichlet_sampler.h"

namespace gelex::detail
{

template <std::size_t ClassCount>
    requires(ClassCount > 1)
class ClassSampler
{
   public:
    struct Posterior
    {
        std::array<double, ClassCount> weights{};
        double total_weight{};
        double log_marginal_kernel{};
    };

    auto begin_sweep(const std::array<double, ClassCount>& probabilities)
        -> void
    {
        uniform_.reset();
        std::ranges::transform(
            probabilities,
            log_probabilities_.begin(),
            [](double probability)
            {
                assert(std::isfinite(probability) && probability > 0.0);
                return std::log(probability);
            });
    }

    [[nodiscard]] auto posterior(
        std::array<double, ClassCount> log_likelihood_kernels) const
        -> Posterior
    {
        for (std::size_t class_index = 0; class_index < ClassCount;
             ++class_index)
        {
            // NOLINTNEXTLINE(cppcoreguidelines-pro-bounds-constant-array-index)
            auto& log_likelihood = log_likelihood_kernels[class_index];
            // NOLINTNEXTLINE(cppcoreguidelines-pro-bounds-constant-array-index)
            log_likelihood += log_probabilities_[class_index];
        }

        const double maximum
            = *std::ranges::max_element(log_likelihood_kernels);
        assert(std::isfinite(maximum));

        double total_weight = 0.0;
        for (double& weight : log_likelihood_kernels)
        {
            weight = std::exp(weight - maximum);
            total_weight += weight;
        }
        assert(std::isfinite(total_weight) && total_weight > 0.0);

        return {
            .weights = log_likelihood_kernels,
            .total_weight = total_weight,
            .log_marginal_kernel = maximum + std::log(total_weight)};
    }

    auto draw(const Posterior& posterior, std::mt19937_64& rng) -> std::size_t
    {
        const double threshold = uniform_(rng) * posterior.total_weight;
        double cumulative_weight = 0.0;
        for (std::size_t class_index = 0; class_index < ClassCount;
             ++class_index)
        {
            // NOLINTNEXTLINE(cppcoreguidelines-pro-bounds-constant-array-index)
            cumulative_weight += posterior.weights[class_index];
            if (threshold < cumulative_weight)
            {
                return class_index;
            }
        }
        return ClassCount - 1;
    }

   private:
    std::array<double, ClassCount> log_probabilities_{};
    std::uniform_real_distribution<double> uniform_{0.0, 1.0};
};

template <MixtureWeightUpdate Update>
class BinaryAllocationUpdater;

template <>
class BinaryAllocationUpdater<MixtureWeightUpdate::Disabled>
{
    using Sampler = ClassSampler<2>;

   public:
    using Posterior = Sampler::Posterior;

    explicit BinaryAllocationUpdater(
        const ProbabilityParameter<
            MixtureWeightUpdate::Disabled>& /*parameter*/) noexcept
    {
    }

    auto begin_sweep(double probability) -> void
    {
        sampler_.begin_sweep({1.0 - probability, probability});
    }

    [[nodiscard]] auto posterior(
        std::array<double, 2> log_likelihood_kernels) const -> Posterior
    {
        return sampler_.posterior(log_likelihood_kernels);
    }

    auto draw(const Posterior& posterior, std::mt19937_64& rng) -> bool
    {
        return sampler_.draw(posterior, rng) == 1;
    }

    auto update(double& /*probability*/, std::mt19937_64& /*rng*/) noexcept
        -> void
    {
    }

   private:
    Sampler sampler_;
};

template <>
class BinaryAllocationUpdater<MixtureWeightUpdate::Enabled>
{
    using Sampler = ClassSampler<2>;

   public:
    using Posterior = Sampler::Posterior;

    explicit BinaryAllocationUpdater(
        const ProbabilityParameter<MixtureWeightUpdate::Enabled>& parameter)
        : beta_{parameter.hyperprior.alpha, parameter.hyperprior.beta}
    {
    }

    auto begin_sweep(double probability) -> void
    {
        sampler_.begin_sweep({1.0 - probability, probability});
        beta_.reset();
        n_success_ = 0;
        n_fail_ = 0;
    }

    [[nodiscard]] auto posterior(
        std::array<double, 2> log_likelihood_kernels) const -> Posterior
    {
        return sampler_.posterior(log_likelihood_kernels);
    }

    auto draw(const Posterior& posterior, std::mt19937_64& rng) -> bool
    {
        const bool success = sampler_.draw(posterior, rng) == 1;
        if (success)
        {
            ++n_success_;
        }
        else
        {
            ++n_fail_;
        }
        return success;
    }

    auto update(double& probability, std::mt19937_64& rng) -> void
    {
        probability = beta_({.n_success = n_success_, .n_fail = n_fail_}, rng);
    }

   private:
    Sampler sampler_;
    BetaSampler<double> beta_;
    Eigen::Index n_success_{};
    Eigen::Index n_fail_{};
};

template <MixtureWeightUpdate Update, std::size_t ClassCount>
class CategoricalAllocationUpdater;

template <std::size_t ClassCount>
class CategoricalAllocationUpdater<MixtureWeightUpdate::Disabled, ClassCount>
{
    using Sampler = ClassSampler<ClassCount>;

   public:
    using Posterior = Sampler::Posterior;

    explicit CategoricalAllocationUpdater(
        const SimplexParameter<
            ClassCount,
            MixtureWeightUpdate::Disabled>& /*parameter*/) noexcept
    {
    }

    auto begin_sweep(const std::array<double, ClassCount>& probabilities)
        -> void
    {
        sampler_.begin_sweep(probabilities);
    }

    [[nodiscard]] auto posterior(
        std::array<double, ClassCount> log_likelihood_kernels) const
        -> Posterior
    {
        return sampler_.posterior(log_likelihood_kernels);
    }

    auto draw(const Posterior& posterior, std::mt19937_64& rng) -> std::size_t
    {
        return sampler_.draw(posterior, rng);
    }

    auto update(
        std::array<double, ClassCount>& /*probabilities*/,
        std::mt19937_64& /*rng*/) noexcept -> void
    {
    }

   private:
    Sampler sampler_;
};

template <std::size_t ClassCount>
class CategoricalAllocationUpdater<MixtureWeightUpdate::Enabled, ClassCount>
{
    static constexpr int class_count = static_cast<int>(ClassCount);
    using Concentrations = Eigen::Vector<double, class_count>;
    using Counts = Eigen::Vector<int, class_count>;
    using Sampler = ClassSampler<ClassCount>;

   public:
    using Posterior = Sampler::Posterior;

    explicit CategoricalAllocationUpdater(
        const SimplexParameter<ClassCount, MixtureWeightUpdate::Enabled>&
            parameter)
        : dirichlet_{Eigen::Map<const Concentrations>{
              parameter.hyperprior.concentration.data()}}
    {
    }

    auto begin_sweep(const std::array<double, ClassCount>& probabilities)
        -> void
    {
        sampler_.begin_sweep(probabilities);
        dirichlet_.reset();
        counts_.setZero();
    }

    [[nodiscard]] auto posterior(
        std::array<double, ClassCount> log_likelihood_kernels) const
        -> Posterior
    {
        return sampler_.posterior(log_likelihood_kernels);
    }

    auto draw(const Posterior& posterior, std::mt19937_64& rng) -> std::size_t
    {
        const std::size_t class_index = sampler_.draw(posterior, rng);
        ++counts_(static_cast<Eigen::Index>(class_index));
        return class_index;
    }

    auto update(
        std::array<double, ClassCount>& probabilities,
        std::mt19937_64& rng) -> void
    {
        const auto sample = dirichlet_(counts_, rng);
        std::copy_n(sample.data(), ClassCount, probabilities.begin());
    }

   private:
    Sampler sampler_;
    DirichletSampler<double> dirichlet_;
    Counts counts_{Counts::Zero()};
};

}  // namespace gelex::detail

#endif  // GELEX_BAYES_GENETIC_KERNEL_ALLOCATION_UPDATER_H_
