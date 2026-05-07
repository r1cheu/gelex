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

#ifndef GELEX_INFRA_LOGGING_PHENO_EVENT_H_
#define GELEX_INFRA_LOGGING_PHENO_EVENT_H_

#include <cstddef>
#include <functional>
#include <optional>
#include <string>
#include <variant>
#include <vector>

namespace gelex
{

struct PhenotypeLoadedEvent
{
    size_t pheno_samples;
    std::string trait_name;
};

struct CovariatesLoadedEvent
{
    std::optional<size_t> num_quantitative_covariates;
    std::optional<size_t> num_discrete_covariates;
    std::vector<std::string> quantitative_names;
    std::vector<std::string> discrete_names;
};

using PhenoEvent = std::variant<PhenotypeLoadedEvent, CovariatesLoadedEvent>;
using PhenoObserver = std::function<void(const PhenoEvent&)>;

}  // namespace gelex

#endif  // GELEX_INFRA_LOGGING_PHENO_EVENT_H_
