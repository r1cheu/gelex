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

#ifndef GELEX_BAYES_GENETIC_KERNEL_DIRICHLET_CONJUGATE_UPDATER_H_
#define GELEX_BAYES_GENETIC_KERNEL_DIRICHLET_CONJUGATE_UPDATER_H_

#include <array>
#include <cstddef>
#include <random>
#include <utility>

#include "gelex/bayes/genetic_family.h"
#include "gelex/bayes/parameter.h"
#include "gelex/bayes/stats/dirichlet_distribution.h"
#include "gelex/bayes/stats/dirichlet_log_kernel.h"

namespace gelex::detail
{

template <std::size_t K, MixtureWeightUpdate Update>
class DirichletConjugateUpdater;

template <std::size_t K>
    requires(K > 1)
class DirichletConjugateUpdater<K, MixtureWeightUpdate::Disabled>
{
   public:
    auto update(
        std::array<double, K>& /*weights*/,
        const std::array<std::size_t, K>& /*counts*/,
        std::mt19937_64& /*rng*/) noexcept -> void
    {
    }
};

template <std::size_t K>
    requires(K > 1)
class DirichletConjugateUpdater<K, MixtureWeightUpdate::Enabled>
{
   public:
    explicit DirichletConjugateUpdater(DirichletLogKernel<K> prior)
        : prior_{std::move(prior)}
    {
    }

    auto update(
        std::array<double, K>& weights,
        const std::array<std::size_t, K>& counts,
        std::mt19937_64& rng) -> void
    {
        const auto posterior = prior_ + make_categorical_likelihood(counts);
        distribution_.reset();
        weights = distribution_(rng, posterior.dirichlet_parameters());
    }

   private:
    DirichletLogKernel<K> prior_;
    DirichletDistribution<K> distribution_;
};

template <std::size_t K, typename T>
    requires(K > 1)
[[nodiscard]] auto make_dirichlet_conjugate_updater(
    const FixedParameter<T>& /*parameter*/)
    -> DirichletConjugateUpdater<K, MixtureWeightUpdate::Disabled>
{
    return {};
}

template <std::size_t K, typename T>
    requires(K > 1)
[[nodiscard]] auto make_dirichlet_conjugate_updater(
    const Parameter<T, DirichletLogKernel<K>>& parameter)
    -> DirichletConjugateUpdater<K, MixtureWeightUpdate::Enabled>
{
    return DirichletConjugateUpdater<K, MixtureWeightUpdate::Enabled>{
        parameter.prior};
}

}  // namespace gelex::detail

#endif  // GELEX_BAYES_GENETIC_KERNEL_DIRICHLET_CONJUGATE_UPDATER_H_
