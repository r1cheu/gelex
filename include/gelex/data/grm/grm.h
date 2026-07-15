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

#ifndef GELEX_DATA_GRM_H_
#define GELEX_DATA_GRM_H_

#include <Eigen/Core>
#include <functional>
#include <span>
#include <string>
#include <string_view>
#include <vector>

#include "gelex/data/genotype_method.h"
#include "gelex/data/marker_range.h"
#include "gelex/infra/logging/grm_event.h"
#include "gelex/types/genetic_mode.h"

namespace gelex
{

class Bed;

struct GrmMatrix
{
    std::string label;  // empty for whole-genome; chromosome name for per-chr
    GeneticMode mode;
    Eigen::MatrixXd grm;
    double denominator;
};

// Builds one GRM per (range, mode), counting each variant once and sharing it
// across modes. Every finished GRM is streamed to a sink rather than held, so
// only the modes of the range in flight are resident.
class GrmBuilder
{
   public:
    using Sink = std::function<void(const GrmMatrix&)>;

    GrmBuilder(
        const Bed& bed,
        GeneticModeSet modes,
        GenotypeMethod method,
        Eigen::Index chunk_size,
        GrmObserver observer = {});

    auto build(std::span<const MarkerRange> ranges, const Sink& sink) -> void;

   private:
    auto accumulate(
        std::string_view label,
        Eigen::Index start,
        Eigen::Index end) -> std::vector<GrmMatrix>;

    const Bed& bed_;
    GeneticModeSet modes_;
    GenotypeMethod method_;
    Eigen::Index chunk_size_;
    GrmObserver observer_;
    Eigen::Index processed_ = 0;
    Eigen::Index total_ = 0;
};

}  // namespace gelex

#endif  // GELEX_DATA_GRM_H_
