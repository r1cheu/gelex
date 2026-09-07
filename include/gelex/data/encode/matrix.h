// Copyright 2026 RuLei Chen
// SPDX-License-Identifier: Apache-2.0

#ifndef GELEX_DATA_ENCODE_MATRIX_H_
#define GELEX_DATA_ENCODE_MATRIX_H_

#include <Eigen/Core>

#include "gelex/data/genotype_method.h"
#include "gelex/genetic_mode.h"

namespace gelex
{

// Encodes an A1-dosage matrix in place; rows are samples, columns are markers,
// and values are 0, 1, 2, or NaN. Invalid loci are replaced with zeros.
auto encode_inplace(
    Eigen::Ref<Eigen::MatrixXd> genotypes,
    GeneticMode mode,
    GenotypeMethod method) -> void;

}  // namespace gelex

#endif  // GELEX_DATA_ENCODE_MATRIX_H_
