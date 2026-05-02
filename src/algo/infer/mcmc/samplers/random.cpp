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

#include "gelex/algo/infer/mcmc/samplers/random.h"

#include <cstddef>

#include <fmt/format.h>
#include <Eigen/Core>

#include "gelex/exception.h"

namespace gelex::mcmc
{

auto RandomSampler::make(const Context& ctx) -> RandomSampler
{
    const auto& effects = ctx.model.random();
    const auto& priors = ctx.priors.random();
    auto& states = ctx.state.random();
    if (effects.size() != priors.size() || effects.size() != states.size())
    {
        throw GelexException(
            fmt::format(
                "random block size mismatch: model={}, priors={}, state={}",
                effects.size(),
                priors.size(),
                states.size()));
    }
    static const bayes::RandomPrior kEmptyPrior{};
    return RandomSampler{Deps{
        .effects = std::span{effects},
        .prior = priors.empty() ? kEmptyPrior : priors.front(),
        .states = std::span{states},
        .residual = ctx.state.residual(),
        .rng = ctx.rng,
    }};
}

auto RandomSampler::sample() -> void
{
    auto& y_adj = deps_.residual.y_adj;
    const double residual_variance = deps_.residual.variance;

    for (std::size_t block = 0; block < deps_.effects.size(); ++block)
    {
        normal_.reset();
        variance_sampler_.reset();

        const auto& effect = deps_.effects[block];
        auto& state = deps_.states[block];

        auto& coeffs = state.coeffs;
        const auto& X = effect.X;
        const auto& XtX_diag = effect.XtX_diag;

        normal_.set_prior_var(state.variance);
        for (Eigen::Index i = 0; i < coeffs.size(); ++i)
        {
            const double old_i = coeffs(i);
            const auto col = X.col(i);
            const double norm = XtX_diag(i);
            const double rhs = col.dot(y_adj) + (norm * old_i);

            const double new_i = normal_(
                stats::NormalSampler<double>::Kernel{
                    .quadratic = norm,
                    .linear = rhs,
                    .scale = residual_variance,
                },
                deps_.rng);
            coeffs(i) = new_i;

            y_adj.array() += (old_i - new_i) * col.array();
        }

        state.variance = variance_sampler_(
            {.n = coeffs.size(), .sum_squares = coeffs.squaredNorm()},
            deps_.rng);
    }
}

}  // namespace gelex::mcmc
