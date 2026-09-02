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

#include "gelex/algo/mcmc/steps/random.h"

#include <Eigen/Core>
#include <cstddef>
#include <fmt/format.h>
#include <random>
#include <span>

#include "gelex/bayes/design.h"
#include "gelex/bayes/prior.h"
#include "gelex/bayes/state.h"
#include "gelex/exception.h"
#include "gelex/infra/stats/normal_sampler.h"

#include "algo/mcmc/invariant.h"

namespace gelex
{

RandomStep::RandomStep(
    const bayes::RandomPrior& prior,
    std::span<const bayes::RandomDesign> designs,
    std::span<bayes::RandomState> states,
    bayes::ResidualState& residual,
    std::mt19937_64& rng)
    : designs_(designs),
      states_(states),
      residual_(residual),
      rng_(rng),
      variance_sampler_(prior.prior())
{
    if (designs.size() != states.size())
    {
        throw GelexException(
            fmt::format(
                "RandomStep: design/state size mismatch: {} != {}",
                designs.size(),
                states.size()));
    }
}

auto RandomStep::step() -> void
{
    const double residual_variance = residual_.variance;
    variance_sampler_.reset();
    for (std::size_t block = 0; block < designs_.size(); ++block)
    {
        const auto& design = designs_[block];
        auto& state = states_[block];
        auto& coeffs = state.coeffs;
        const auto& X = design.X;
        const auto& XtX_diag = design.XtX_diag;

        normal_.reset();
        normal_.set_prior_var(state.variance);
        for (Eigen::Index i = 0; i < coeffs.size(); ++i)
        {
            const auto col = X.col(i);
            const double old_i = coeffs(i);
            DenseResidualAdjustmentGuard guard{col, coeffs(i), residual_};
            const double rhs = col.dot(residual_.y_adj) + (XtX_diag(i) * old_i);
            coeffs(i) = normal_(
                NormalSampler<double>::Kernel{
                    .quadratic = XtX_diag(i),
                    .linear = rhs,
                    .scale = residual_variance,
                },
                rng_);
        }
        state.variance = variance_sampler_(
            {.n = coeffs.size(), .sum_squares = coeffs.squaredNorm()}, rng_);
    }
}

}  // namespace gelex
