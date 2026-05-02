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

#ifndef GELEX_MODEL_BAYES_GENOTYPE_STORAGE_H_
#define GELEX_MODEL_BAYES_GENOTYPE_STORAGE_H_

#include <variant>

#include <Eigen/Core>

#include "gelex/data/genotype/genotype_matrix.h"
#include "gelex/data/genotype/genotype_mmap.h"

namespace gelex::bayes
{

using GenotypeStorage
    = std::variant<genotype::GenotypeMap, genotype::GenotypeMatrix>;

inline auto get_matrix_ref(const GenotypeStorage& storage)
    -> Eigen::Ref<const Eigen::MatrixXd>
{
    return std::visit(
        [](const auto& s) -> Eigen::Ref<const Eigen::MatrixXd>
        { return s.matrix(); },
        storage);
}

inline auto get_rows(const GenotypeStorage& storage) -> Eigen::Index
{
    return std::visit([](const auto& s) { return s.rows(); }, storage);
}

inline auto get_cols(const GenotypeStorage& storage) -> Eigen::Index
{
    return std::visit([](const auto& s) { return s.cols(); }, storage);
}

inline auto get_means(const GenotypeStorage& storage) -> const Eigen::VectorXd&
{
    return std::visit(
        [](const auto& s) -> const Eigen::VectorXd& { return s.mean(); },
        storage);
}

inline auto get_stddev(const GenotypeStorage& storage) -> const Eigen::VectorXd&
{
    return std::visit(
        [](const auto& s) -> const Eigen::VectorXd& { return s.stddev(); },
        storage);
}

inline auto is_monomorphic_variant(
    const GenotypeStorage& storage,
    Eigen::Index idx) -> bool
{
    return std::visit(
        [idx](const auto& s) { return s.is_monomorphic(idx); }, storage);
}

inline auto num_mono_variant(const GenotypeStorage& storage) -> Eigen::Index
{
    return std::visit([](const auto& s) { return s.num_mono(); }, storage);
}

}  // namespace gelex::bayes

#endif  // GELEX_MODEL_BAYES_GENOTYPE_STORAGE_H_
