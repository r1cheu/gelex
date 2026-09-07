// Copyright 2026 RuLei Chen
// SPDX-License-Identifier: Apache-2.0

#ifndef GELEX_DATA_SNP_LUT_H_
#define GELEX_DATA_SNP_LUT_H_

#include <Eigen/Core>

namespace gelex
{

// Indexed by raw BED code: 00=A1A1, 01=missing, 10=A1A2, 11=A2A2.
using SnpLut = Eigen::Array4d;
using SnpLutMatrix = Eigen::Array<double, 4, Eigen::Dynamic, Eigen::ColMajor>;

}  // namespace gelex

#endif  // GELEX_DATA_SNP_LUT_H_
