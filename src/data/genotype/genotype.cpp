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

#include "gelex/data/genotype/genotype.h"

#include <algorithm>
#include <cstdint>
#include <type_traits>
#include <variant>

#include <Eigen/Core>
#include <vector>

namespace gelex::genotype
{

auto Genotype::matrix() const noexcept -> Eigen::Ref<const Eigen::MatrixXd>
{
    return std::visit(
        [](const auto& s) -> Eigen::Ref<const Eigen::MatrixXd>
        {
            using S = std::decay_t<decltype(s)>;
            if constexpr (std::is_same_v<S, OwnedStorage>)
            {
                return s.data;
            }
            else
            {
                return s.view;
            }
        },
        storage_);
}

auto Genotype::mean() const noexcept -> const Eigen::VectorXd&
{
    return std::visit(
        [](const auto& s) -> const Eigen::VectorXd& { return s.mean; },
        storage_);
}

auto Genotype::stddev() const noexcept -> const Eigen::VectorXd&
{
    return std::visit(
        [](const auto& s) -> const Eigen::VectorXd& { return s.stddev; },
        storage_);
}

auto Genotype::allele_freq() const noexcept -> const Eigen::VectorXd&
{
    return std::visit(
        [](const auto& s) -> const Eigen::VectorXd& { return s.allele_freq; },
        storage_);
}

auto Genotype::mono_indices() const noexcept -> const std::vector<int64_t>&
{
    return std::visit(
        [](const auto& s) -> const std::vector<int64_t>&
        { return s.mono_indices; },
        storage_);
}

auto Genotype::num_mono() const noexcept -> int64_t
{
    return static_cast<int64_t>(mono_indices().size());
}

auto Genotype::rows() const noexcept -> int64_t
{
    return static_cast<int64_t>(matrix().rows());
}

auto Genotype::cols() const noexcept -> int64_t
{
    return static_cast<int64_t>(matrix().cols());
}

auto Genotype::is_monomorphic(Eigen::Index marker_idx) const noexcept -> bool
{
    const auto& v = mono_indices();
    return std::ranges::binary_search(v, static_cast<int64_t>(marker_idx));
}

}  // namespace gelex::genotype
