// Copyright 2026 RuLei Chen
// SPDX-License-Identifier: Apache-2.0

#ifndef GELEX_FREQ_REML_IO_H_
#define GELEX_FREQ_REML_IO_H_

#include <span>
#include <string>
#include <string_view>
#include <vector>

#include "gelex/freq/reml/summary.h"

namespace gelex
{

class FreqModel;
class FreqState;

auto write_summary(
    const FreqModel& model,
    const FreqState& state,
    double loglike,
    std::string_view prefix) -> void;

auto write_effects(
    const FreqModel& model,
    const FreqState& state,
    std::span<const std::string> sample_ids,
    std::string_view prefix) -> void;

auto write_loco_summary(
    const std::vector<LocoRemlResult>& results,
    std::string_view prefix) -> void;

}  // namespace gelex

#endif  // GELEX_FREQ_REML_IO_H_
