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

#ifndef GELEX_INFRA_LOGGING_ASSOC_EVENT_H_
#define GELEX_INFRA_LOGGING_ASSOC_EVENT_H_

#include <cstddef>
#include <functional>
#include <string>
#include <variant>
#include <vector>

#include "gelex/data/genotype/genotype_processor.h"
#include "gelex/infra/logging/reml_event.h"
#include "gelex/types/genetic_effect_type.h"

namespace gelex
{

struct AssocConfigLoadedEvent
{
    ModelType model_type;
    bool loco;

    GenotypeProcessMethod geno_method;

    int max_iter;
    double tol;
};

struct AssocRemlStartedEvent
{
    std::string chr_name;
};

struct AssocScanSummaryEvent
{
    size_t total_snps;
    int chunk_size;
    bool loco;
};

struct AssocScanProgressEvent
{
    size_t current;
    size_t total;
};

struct AssocLocoPhaseEvent
{
    std::string chr_name;
    std::string phase;
};

struct AssocLocoRemlSummaryEvent
{
    std::vector<LocoRemlResult> results;
};

struct AssocCompleteEvent
{
    std::string out_prefix;
};

using AssocEvent = std::variant<
    AssocConfigLoadedEvent,
    AssocRemlStartedEvent,
    AssocScanSummaryEvent,
    AssocScanProgressEvent,
    AssocLocoPhaseEvent,
    AssocLocoRemlSummaryEvent,
    AssocCompleteEvent>;

using AssocObserver = std::function<void(const AssocEvent&)>;

}  // namespace gelex

#endif  // GELEX_INFRA_LOGGING_ASSOC_EVENT_H_
