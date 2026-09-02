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
