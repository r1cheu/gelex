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

#ifndef GELEX_INFRA_LOGGING_PREDICT_EVENT_H_
#define GELEX_INFRA_LOGGING_PREDICT_EVENT_H_

#include <cstddef>
#include <functional>
#include <string>
#include <variant>

#include "gelex/data/genotype/process_method.h"

namespace gelex
{

struct PredictBannerEvent
{
};

struct PredictParamsLoadedEvent
{
    std::string_view bfile_prefix;
    std::string_view gfile_prefix;
    GenotypeProcessMethod geno_method{};
};

struct PredictSnpSelectionEvent
{
    size_t num_matched{};
    size_t num_missing{};
    size_t num_mismatched{};
    size_t num_total{};
    std::string bfile_path;
    std::string snp_effect_path;
};

struct PredictDataLoadedEvent
{
    size_t num_samples{};
    size_t num_snps{};
    size_t num_covar_terms{};
};

struct PredictResultsWrittenEvent
{
    std::string output_path;
    size_t num_samples{};
};

using PredictEvent = std::variant<
    PredictBannerEvent,
    PredictParamsLoadedEvent,
    PredictSnpSelectionEvent,
    PredictDataLoadedEvent,
    PredictResultsWrittenEvent>;
using PredictObserver = std::function<void(const PredictEvent&)>;

}  // namespace gelex

#endif  // GELEX_INFRA_LOGGING_PREDICT_EVENT_H_
