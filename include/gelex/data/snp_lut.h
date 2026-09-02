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
