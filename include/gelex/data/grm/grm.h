// Copyright 2026 RuLei Chen
// SPDX-License-Identifier: Apache-2.0

#ifndef GELEX_DATA_GRM_GRM_H_
#define GELEX_DATA_GRM_GRM_H_

#include <Eigen/Core>
#include <cstddef>
#include <functional>
#include <span>
#include <string>
#include <string_view>
#include <vector>

#include "gelex/data/genotype_method.h"
#include "gelex/data/marker_range.h"
#include "gelex/genetic_mode.h"

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
        std::function<void(std::size_t)> observer = {});

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
    std::function<void(std::size_t)> observer_;
    Eigen::Index processed_ = 0;
};

}  // namespace gelex

#endif  // GELEX_DATA_GRM_GRM_H_
