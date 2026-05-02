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

#include "gelex/types/vi_result.h"

#include "gelex/model/bayes/genotype_storage.h"
#include "gelex/model/bayes/model.h"

namespace gelex
{

vi::Result::Result(vi::State&& state, const BayesModel& model)
    : fixed_(FixedSummary(FixedSamples(model.fixed())))
{
    const double res_var = state.residual().variance;

    auto& fc = fixed_.coeffs;
    fc.mean = std::move(state.fixed().coeffs);
    fc.stddev = (res_var / model.fixed().XtX_diag.array()).sqrt();

    const auto& random_effects = model.random();
    for (std::size_t i = 0; i < random_effects.size(); ++i)
    {
        RandomSummary rs{RandomSamples(random_effects[i])};
        auto& rs_state = state.random()[i];

        rs.coeffs.mean = std::move(rs_state.coeffs);
        const auto inv_scaler = 1.0
                                / (random_effects[i].XtX_diag.array()
                                   + res_var / rs_state.variance);
        rs.coeffs.stddev = (res_var * inv_scaler).sqrt();

        rs.variance.mean = Eigen::VectorXd::Constant(1, rs_state.variance);
        rs.variance.stddev = Eigen::VectorXd::Zero(1);

        random_.push_back(std::move(rs));
    }

    for (auto& gs : state.genetics())
    {
        genetics_.push_back(
            vi::GeneticSummary{
                .type = gs.type,
                .coeffs = std::move(gs.coeffs),
            });
    }

    if (const auto* effect = model.genetic(GeneticMode::A); effect)
    {
        p_freq_ = bayes::get_means(effect->X).array() / 2;
    }
}

}  // namespace gelex
