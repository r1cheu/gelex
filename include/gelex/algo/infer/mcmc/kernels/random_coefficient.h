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

#ifndef GELEX_ALGO_INFER_MCMC_KERNELS_RANDOM_COEFFICIENT_H_
#define GELEX_ALGO_INFER_MCMC_KERNELS_RANDOM_COEFFICIENT_H_

#include <random>

#include "gelex/infra/stats/conjugate_prior.h"

namespace gelex::mcmc
{

class RandomCoefficientKernel
{
   public:
    auto prepare(double prior_variance) -> void
    {
        normal_.reset();
        normal_.set_prior_var(prior_variance);
    }

    auto sample(
        double xtx_diag_i,
        double rhs,
        double residual_variance,
        std::mt19937_64& rng) -> double
    {
        return normal_(
            stats::NormalSampler<double>::Kernel{
                .quadratic = xtx_diag_i,
                .linear = rhs,
                .scale = residual_variance,
            },
            rng);
    }

   private:
    stats::NormalSampler<double> normal_{0.0};
};

}  // namespace gelex::mcmc

#endif  // GELEX_ALGO_INFER_MCMC_KERNELS_RANDOM_COEFFICIENT_H_
