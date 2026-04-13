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

#include "gelex/pipeline/vi_engine.h"

#include <utility>

#include <fmt/format.h>

#include "gelex/algo/infer/vi/solver.h"
#include "gelex/algo/infer/vi/trait_model.h"
#include "gelex/exception.h"
#include "gelex/infra/logging/notify.h"
#include "gelex/model/bayes/model.h"
#include "gelex/model/bayes/writer/vi_result_writer.h"
#include "gelex/pipeline/model_builder.h"

namespace gelex
{

vi::FitEngine::FitEngine(Config config) : config_(std::move(config)) {}

auto vi::FitEngine::run(
    PhenoPipe&& pheno,
    GenoPipe&& geno,
    const FitObserver& observer) -> void
{
    if (config_.method.base != BayesBase::RR)
    {
        throw GelexException(
            fmt::format("CAVI only supports BayesRR, got: {}", config_.method));
    }

    auto model = build_bayes_model(std::move(pheno), std::move(geno));
    auto priors = build_bayes_priors(
        PriorOverrides{
            .method = config_.method,
            .phenotype_variance = model.phenotype_variance(),
            .pi = std::nullopt,
            .dpi = std::nullopt,
            .multiplier = std::nullopt,
            .dmultiplier = std::nullopt,
            .positive_prob = 0.5,
        },
        model);

    vi::Result result = [&]
    {
        if (config_.method.dominance)
        {
            vi::Solver engine(config_.params, vi::RRd{});
            return engine.run(model, priors, observer);
        }
        vi::Solver engine(config_.params, vi::RR{});
        return engine.run(model, priors, observer);
    }();

    vi::ResultWriter writer(result, config_.bfile_prefix + ".bim");
    writer.save(config_.out_prefix);

    notify(observer, FitResultsSavedEvent{.out_prefix = config_.out_prefix});
}

}  // namespace gelex
