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

#ifndef GELEX_TYPES_SNP_EFFECTS_H
#define GELEX_TYPES_SNP_EFFECTS_H

#include <cassert>
#include <optional>
#include <vector>

#include <Eigen/Core>

#include "gelex/types/snp_index.h"

namespace gelex
{

class SnpEffects
{
   public:
    explicit SnpEffects(size_t initial_capacity = 0);

    auto emplace(SnpMeta meta, double additive_effect, double A1_frequency)
        -> void
    {
        assert(dominance_data_.empty() && "Error: Mixed usage! Dominance mode is active but 2-arg method called.");

        index_.emplace(std::move(meta));
        additive_data_.push_back(additive_effect);
        frequencies_data_.push_back(A1_frequency);
    }

    auto emplace(
        SnpMeta meta,
        double additive_effect,
        double dominance_effect,
        double A1_frequency) -> void
    {
        assert((additive_data_.empty() || !dominance_data_.empty()) &&
               "Error: Mixed usage! Additive-only mode is active but 4-arg method called.");

        index_.emplace(std::move(meta));
        additive_data_.push_back(additive_effect);
        dominance_data_.push_back(dominance_effect);
        frequencies_data_.push_back(A1_frequency);
    }

    auto index() -> SnpIndex& { return index_; }
    auto index() const -> const SnpIndex& { return index_; }

    auto operator[](size_t i) -> SnpMeta& { return index_[i]; }
    auto operator[](size_t i) const -> const SnpMeta& { return index_[i]; }

    auto find_index(std::string_view snp_id) const
        -> std::optional<Eigen::Index>
    {
        return index_.find_index(snp_id);
    }

    [[nodiscard]] auto size() const -> size_t { return index_.size(); }

    auto additive_effects() -> Eigen::Map<Eigen::VectorXd>
    {
        return {
            additive_data_.data(),
            static_cast<Eigen::Index>(additive_data_.size())};
    }
    auto additive_effects() const -> Eigen::Map<const Eigen::VectorXd>
    {
        return {
            additive_data_.data(),
            static_cast<Eigen::Index>(additive_data_.size())};
    }

    auto dominance_effects() -> Eigen::Map<Eigen::VectorXd>
    {
        return {
            dominance_data_.data(),
            static_cast<Eigen::Index>(dominance_data_.size())};
    }
    auto dominance_effects() const -> Eigen::Map<const Eigen::VectorXd>
    {
        return {
            dominance_data_.data(),
            static_cast<Eigen::Index>(dominance_data_.size())};
    }

    auto frequencies() -> Eigen::Map<Eigen::VectorXd>
    {
        return {
            frequencies_data_.data(),
            static_cast<Eigen::Index>(frequencies_data_.size())};
    }
    auto frequencies() const -> Eigen::Map<const Eigen::VectorXd>
    {
        return {
            frequencies_data_.data(),
            static_cast<Eigen::Index>(frequencies_data_.size())};
    }

    auto shrink_to_fit() -> void;
    auto clear() -> void;

   private:
    SnpIndex index_;

    std::vector<double> additive_data_;
    std::vector<double> dominance_data_;
    std::vector<double> frequencies_data_;
};

}  // namespace gelex
#endif  // GELEX_TYPES_SNP_EFFECTS_H
