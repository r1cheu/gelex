// Copyright 2026 RuLei Chen
// SPDX-License-Identifier: Apache-2.0

#ifndef GELEX_BAYES_DETAIL_NORMAL_VARIANCE_CONJUGATE_UPDATER_H_
#define GELEX_BAYES_DETAIL_NORMAL_VARIANCE_CONJUGATE_UPDATER_H_

#include <cstddef>
#include <random>

#include "gelex/bayes/stats/scaled_inv_chi2_distribution.h"
#include "gelex/bayes/stats/scaled_inv_chi2_log_kernel.h"

namespace gelex::detail
{

class NormalVarianceConjugateUpdater
{
   public:
    explicit NormalVarianceConjugateUpdater(ScaledInvChi2LogKernel prior)
        : prior_{prior}
    {
    }

    auto update(
        double& variance,
        std::size_t count,
        double sum_squares,
        std::mt19937_64& rng) -> void
    {
        const auto likelihood
            = make_normal_variance_likelihood(count, sum_squares);
        const auto posterior = prior_ + likelihood;
        distribution_.reset();
        variance = distribution_(rng, posterior.scaled_inv_chi2_parameters());
    }

   private:
    ScaledInvChi2LogKernel prior_;
    ScaledInvChi2Distribution<> distribution_;
};

}  // namespace gelex::detail

#endif  // GELEX_BAYES_DETAIL_NORMAL_VARIANCE_CONJUGATE_UPDATER_H_
