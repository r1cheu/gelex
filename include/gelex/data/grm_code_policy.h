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

#ifndef GELEX_DATA_GRM_CODE_POLICY_H
#define GELEX_DATA_GRM_CODE_POLICY_H

#include <Eigen/Core>

namespace gelex::grm
{

constexpr double EPSILON = 1e-10;

namespace detail
{
auto additive_mean_center(
    Eigen::Ref<Eigen::MatrixXd> genotype,
    Eigen::VectorXd* freqs = nullptr) -> void;
}  // namespace detail

struct Su
{
    auto operator()(
        Eigen::Ref<Eigen::MatrixXd> genotype,
        bool use_additive,
        Eigen::VectorXd* freqs = nullptr) const -> void;
};

struct Zeng
{
    auto operator()(
        Eigen::Ref<Eigen::MatrixXd> genotype,
        bool use_additive,
        Eigen::VectorXd* freqs = nullptr) const -> void;
};

struct Yang
{
    auto operator()(
        Eigen::Ref<Eigen::MatrixXd> genotype,
        bool use_additive,
        Eigen::VectorXd* freqs = nullptr) const -> void;
};

struct Vitezica
{
    auto operator()(
        Eigen::Ref<Eigen::MatrixXd> genotype,
        bool use_additive,
        Eigen::VectorXd* freqs = nullptr) const -> void;
};

}  // namespace gelex::grm

#endif  // GELEX_DATA_GRM_CODE_POLICY_H
