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

#ifndef GELEX_DATA_SNP_ALIGNMENT_H_
#define GELEX_DATA_SNP_ALIGNMENT_H_

#include <Eigen/Core>
#include <string>
#include <vector>

#include "gelex/data/dataframe/dataframe.h"

namespace gelex
{

class Bed;
struct SnpStats;

inline constexpr double MAX_SNP_MISSING_RATIO = 0.2;

// Maps training SNP effects (canonical axis) onto the columns of a prediction
// .bim by SNP id, orienting each match to the training A1. Three outcomes per
// training SNP: same orientation, allele-swapped (encoding swaps
// code[0]/code[2], equivalent to 2 - x), or missing (id absent, or alleles
// neither match nor swap). Pure over allele tables: no genotype I/O, no
// imputation.
struct AlignmentPlan
{
    std::vector<Eigen::Index> source_col;  // bim columns to read (dense)
    std::vector<Eigen::Index> train_pos;   // parallel: training-length target
    std::vector<char> flip;  // parallel: swap code[0]/code[2] when nonzero
    std::vector<Eigen::Index> missing_pos;  // training positions with no source

    Eigen::Index train_count{};  // width of the training axis (output columns)
    Eigen::Index num_same{};
    Eigen::Index num_flip{};
    Eigen::Index num_absent{};  // id not present in bim
    Eigen::Index
        num_incompatible{};  // id present but alleles neither match nor swap
};

[[nodiscard]] auto build_snp_alignment(
    const DataFrame<std::string>& snp_effects,
    const DataFrame<std::string>& bim_df) -> AlignmentPlan;

// Executes an AlignmentPlan against opened genotypes for one training mode:
// expands each planned bim column straight from its packed form using the
// training codes (stats.code), swapping code[0]/code[2] for flipped SNPs, and
// scatters onto a plan.train_count-wide matrix in training order. Columns for
// missing_pos are filled with the encoding's missing value (zero contribution).
[[nodiscard]] auto expand_aligned_genotypes(
    const Bed& bed,
    const AlignmentPlan& plan,
    const SnpStats& stats) -> Eigen::MatrixXd;

}  // namespace gelex

#endif  // GELEX_DATA_SNP_ALIGNMENT_H_
