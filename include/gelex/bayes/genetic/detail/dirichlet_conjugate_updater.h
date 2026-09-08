// Copyright 2026 RuLei Chen
// SPDX-License-Identifier: Apache-2.0

#ifndef GELEX_BAYES_GENETIC_DETAIL_DIRICHLET_CONJUGATE_UPDATER_H_
#define GELEX_BAYES_GENETIC_DETAIL_DIRICHLET_CONJUGATE_UPDATER_H_

#include <array>
#include <cstddef>
#include <random>
#include <utility>

#include "gelex/bayes/stats/dirichlet_distribution.h"
#include "gelex/bayes/stats/dirichlet_log_kernel.h"

namespace gelex::detail
{

template <std::size_t K>
    requires(K > 1)
class DirichletConjugateUpdater
{
   public:
    explicit DirichletConjugateUpdater(DirichletLogKernel<K> prior)
        : prior_{std::move(prior)}
    {
    }

    [[nodiscard]] auto draw(
        const std::array<std::size_t, K>& counts,
        std::mt19937_64& rng) -> std::array<double, K>
    {
        const auto posterior = prior_ + make_categorical_likelihood(counts);
        distribution_.reset();
        return distribution_(rng, posterior.dirichlet_parameters());
    }

   private:
    DirichletLogKernel<K> prior_;
    DirichletDistribution<K> distribution_;
};

}  // namespace gelex::detail

#endif  // GELEX_BAYES_GENETIC_DETAIL_DIRICHLET_CONJUGATE_UPDATER_H_
