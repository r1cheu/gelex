// Copyright 2026 RuLei Chen
// SPDX-License-Identifier: Apache-2.0

#ifndef GELEX_FREQ_REML_PROGRESS_H_
#define GELEX_FREQ_REML_PROGRESS_H_

#include <cstddef>
#include <functional>
#include <string>
#include <variant>
#include <vector>

namespace gelex
{

struct RemlIterationEvent
{
    size_t iter;
    double loglike;
    std::vector<std::string> labels;
    std::vector<double> variances;
};

struct RemlConstrainedEvent
{
    size_t num_constrained;
    size_t num_total;
};

using RemlEvent = std::variant<RemlIterationEvent, RemlConstrainedEvent>;
using RemlObserver = std::function<void(const RemlEvent&)>;

}  // namespace gelex

#endif  // GELEX_FREQ_REML_PROGRESS_H_
