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

#ifndef GELEX_ALGO_INFER_MCMC_KERNELS_CONCEPT_H_
#define GELEX_ALGO_INFER_MCMC_KERNELS_CONCEPT_H_

#include <concepts>
#include <random>
#include <utility>

#include <Eigen/Core>

namespace gelex::mcmc
{

template <typename K>
concept GeneticKernel = requires(
    K& k,
    Eigen::Index marker_index,
    double xtx_diag,
    double rhs,
    double residual_variance,
    std::mt19937_64& rng) {
    { k.prepare() } -> std::same_as<void>;
    {
        k.sample(marker_index, xtx_diag, rhs, residual_variance, rng)
    } -> std::same_as<double>;
    { k.commit(rng) } -> std::same_as<void>;
};

template <typename K>
concept GeneticJointKernel = requires(
    K& k,
    Eigen::Index marker_index,
    double first_xtx_diag,
    double first_rhs,
    double second_xtx_diag,
    double second_rhs,
    double residual_variance,
    std::mt19937_64& rng) {
    { k.prepare() } -> std::same_as<void>;
    {
        k.sample(
            marker_index,
            first_xtx_diag,
            first_rhs,
            second_xtx_diag,
            second_rhs,
            residual_variance,
            rng)
    } -> std::same_as<std::pair<double, double>>;
    { k.commit(rng) } -> std::same_as<void>;
};

}  // namespace gelex::mcmc

#endif  // GELEX_ALGO_INFER_MCMC_KERNELS_CONCEPT_H_
