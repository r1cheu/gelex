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

#ifndef GELEX_ALGO_INFER_MCMC_KERNELS_RESIDUAL_VARIANCE_H_
#define GELEX_ALGO_INFER_MCMC_KERNELS_RESIDUAL_VARIANCE_H_

#include <random>

#include <Eigen/Core>

#include "gelex/infra/stats/conjugate_prior.h"
#include "gelex/model/bayes/prior.h"

namespace gelex::mcmc
{

class ResidualVarianceKernel
{
   public:
    explicit ResidualVarianceKernel(const bayes::ResidualPrior& prior)
        : sampler_(prior.prior().degrees_of_freedom(), prior.prior().scale())
    {
    }

    auto prepare() -> void { sampler_.reset(); }

    auto sample(
        Eigen::Index num_individuals,
        const Eigen::VectorXd& y_adj,
        std::mt19937_64& rng) -> double
    {
        return sampler_(
            {.n = num_individuals, .sum_squares = y_adj.squaredNorm()}, rng);
    }

   private:
    stats::ScaledInvChi2Sampler<double> sampler_;
};

}  // namespace gelex::mcmc

#endif  // GELEX_ALGO_INFER_MCMC_KERNELS_RESIDUAL_VARIANCE_H_
