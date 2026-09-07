// Copyright 2026 RuLei Chen
// SPDX-License-Identifier: Apache-2.0

#ifndef GELEX_FREQ_REML_SUMMARY_H_
#define GELEX_FREQ_REML_SUMMARY_H_

#include <cstddef>
#include <string>
#include <vector>

#include "gelex/freq/reml/operators.h"

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

#endif  // GELEX_FREQ_REML_SUMMARY_H_
