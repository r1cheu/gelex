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

#ifndef GELEX_ALGO_INFER_VI_KERNELS_CONCEPT_H_
#define GELEX_ALGO_INFER_VI_KERNELS_CONCEPT_H_

#include <concepts>

#include <Eigen/Core>

namespace gelex::vi
{

struct CaviUpdate
{
    double mean{};
    double sigma2{};
};

template <typename K>
concept GeneticKernel = requires(
    K& k,
    Eigen::Index marker_index,
    double xtx_diag,
    double rhs,
    double residual_variance) {
    { k.prepare() } -> std::same_as<void>;
    {
        k.update(marker_index, xtx_diag, rhs, residual_variance)
    } -> std::same_as<CaviUpdate>;
    { k.commit() } -> std::same_as<void>;
};

}  // namespace gelex::vi

#endif  // GELEX_ALGO_INFER_VI_KERNELS_CONCEPT_H_
