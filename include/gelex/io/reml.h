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

#ifndef GELEX_IO_REML_H_
#define GELEX_IO_REML_H_

#include <span>
#include <string>
#include <string_view>
#include <vector>

#include "gelex/algo/reml/loco_result.h"

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

#endif  // GELEX_IO_REML_H_
