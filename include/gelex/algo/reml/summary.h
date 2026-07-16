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

#ifndef GELEX_ALGO_REML_SUMMARY_H_
#define GELEX_ALGO_REML_SUMMARY_H_

#include <cstddef>
#include <string>
#include <vector>

#include "gelex/algo/reml/operators.h"

namespace gelex
{

struct VarianceComponent
{
    std::string name;
    double variance{};
    double variance_se{};
    double variance_ratio{};
    double variance_ratio_se{};
    bool at_boundary{};  // clamped to the constraint floor; Wald test invalid
};

// Reportable scalar outcome of a single REML fit. Copyable and decoupled from
// the live FreqState so callers (e.g. per-chromosome LOCO) can snapshot it.
struct RemlSummary
{
    double loglike{};
    bool converged{};
    std::size_t iter_count{};
    std::vector<VarianceComponent> random;
    double residual_variance{};
    double residual_variance_se{};
};

// Self-contained outcome of one fit: the reportable scalar summary plus the
// heavy GWAS projection operators. Callers take whichever part they need.
struct RemlFit
{
    RemlSummary summary;
    GwasOperators operators;
};

// One chromosome's REML fit summary in a LOCO scan: the leave-one-chromosome
// -out variance components tagged with the held-out chromosome.
struct LocoRemlResult
{
    std::string chr_name;
    RemlSummary summary;
};

}  // namespace gelex

#endif  // GELEX_ALGO_REML_SUMMARY_H_
