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

#include <string>
#include <vector>

#include <Eigen/Core>

#include "gelex/data/dataframe/dataframe.h"

namespace gelex
{

class Bed;

inline constexpr double MAX_SNP_MISSING_RATIO = 0.2;

// Maps training SNP effects (canonical axis) onto the columns of a prediction
// .bim by SNP id, orienting each match to the training A1. Three outcomes per
// training SNP: same orientation, allele-swapped (dosage needs 2 - x), or
// missing (id absent, or alleles neither match nor swap). Pure over allele
// tables: no genotype I/O, no imputation.
struct AlignmentPlan
{
    std::vector<Eigen::Index> source_col;  // bim columns to read (dense)
    std::vector<Eigen::Index> train_pos;   // parallel: training-length target
    std::vector<char> flip;                // parallel: apply 2 - x when nonzero
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

// Executes an AlignmentPlan against opened genotypes: reads the planned bim
// columns, orients flipped SNPs (2 - x), and scatters them onto a
// plan.train_count-wide matrix in training order. Columns for missing_pos stay
// NaN, which the encoding lookup later maps to a zero contribution.
[[nodiscard]] auto load_aligned_genotypes(
    const Bed& bed,
    const AlignmentPlan& plan) -> Eigen::MatrixXd;

}  // namespace gelex

#endif  // GELEX_DATA_SNP_ALIGNMENT_H_
