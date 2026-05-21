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

#include "gelex/engine/vi.h"

#include <utility>

#include <fmt/format.h>

#include "gelex/algo/infer/vi/recipes.h"
#include "gelex/algo/infer/vi/solver.h"
#include "gelex/exception.h"
#include "gelex/infra/logging/notify.h"
#include "gelex/io/vi/result_writer.h"
#include "gelex/model/bayes/legacy_algorithm_shape.h"
#include "gelex/model/bayes/legacy_bayes_policy.h"
#include "gelex/model/bayes/model.h"

namespace gelex
{

vi::Engine::Engine(Config config) : config_(std::move(config)) {}

auto vi::Engine::run(
    const BayesModel& model,
    bayes::LegacyBayesMethod method,
    const VIObserver& observer) -> void
{
    if (config_.method.base != BayesBase::RR)
    {
        throw GelexException(
            fmt::format("CAVI only supports BayesRR, got: {}", config_.method));
    }

    const auto shape = bayes::resolve_shape(
        bayes::policy_for(config_.method.base), config_.requested_effects);

    vi::Result result = [&]
    {
        switch (shape)
        {
            case bayes::AlgorithmShape::a_only:
            {
                vi::Solver engine(
                    config_.params,
                    vi::make_bayes_rr_chain<bayes::AlgorithmShape::a_only>);
                return engine.run(model, method, observer);
            }
            case bayes::AlgorithmShape::d_only:
            {
                vi::Solver engine(
                    config_.params,
                    vi::make_bayes_rr_chain<bayes::AlgorithmShape::d_only>);
                return engine.run(model, method, observer);
            }
            case bayes::AlgorithmShape::ad_independent:
            {
                vi::Solver engine(
                    config_.params,
                    vi::make_bayes_rr_chain<
                        bayes::AlgorithmShape::ad_independent>);
                return engine.run(model, method, observer);
            }
            case bayes::AlgorithmShape::ad_joint:
                throw GelexException("BayesRR does not support ad_joint shape");
        }
        std::unreachable();
    }();

    vi::ResultWriter writer(result, config_.bfile_prefix + ".bim");
    writer.save(config_.out_prefix);

    notify(observer, FitResultsSavedEvent{.out_prefix = config_.out_prefix});
}

}  // namespace gelex
