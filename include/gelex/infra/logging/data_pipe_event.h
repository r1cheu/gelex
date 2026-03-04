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

#ifndef GELEX_INFRA_LOGGING_DATA_PIPE_EVENT_H_
#define GELEX_INFRA_LOGGING_DATA_PIPE_EVENT_H_

#include <functional>
#include <string>
#include <variant>
#include <vector>

#include "gelex/types/freq_effect.h"

namespace gelex
{

struct DataPipeSectionEvent
{
};

struct PhenotypeLoadedEvent
{
    size_t pheno_samples;
    std::string trait_name;
    size_t geno_samples;
};

struct CovariatesLoadedEvent
{
    std::optional<size_t> num_quantitative_covariates;
    std::optional<size_t> num_discrete_covariates;
    std::vector<std::string> quantitative_names;
    std::vector<std::string> discrete_names;
};

struct IntersectionEvent
{
    size_t common_samples;
    size_t excluded_samples;
};

struct GenotypeLoadedEvent
{
    bool is_dominance;  // false = Additive, true = Dominance
    int64_t num_snps;
    int64_t monomorphic_snps;
};

struct GrmLoadedEvent
{
    size_t num_samples;
    freq::GrmType type;
};

using DataPipeEvent = std::variant<
    DataPipeSectionEvent,
    PhenotypeLoadedEvent,
    CovariatesLoadedEvent,
    IntersectionEvent,
    GenotypeLoadedEvent,
    GrmLoadedEvent>;

using DataPipeObserver = std::function<void(const DataPipeEvent&)>;

}  // namespace gelex

#endif  // GELEX_INFRA_LOGGING_DATA_PIPE_EVENT_H_
