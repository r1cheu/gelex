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

#ifndef GELEX_TYPES_SNP_INFO_H
#define GELEX_TYPES_SNP_INFO_H

#include <cassert>
#include <optional>
#include <string>
#include <string_view>
#include <unordered_map>
#include <vector>

#include <Eigen/Core>

namespace gelex
{

struct SnpMeta
{
    std::string chrom;
    std::string id;
    int pos;
    char A1;
    char A2;
};

class SnpEffects
{
   public:
    using VectorMap = Eigen::Map<Eigen::VectorXd>;
    using ConstVectorMap = Eigen::Map<const Eigen::VectorXd>;

    using iterator = std::vector<SnpMeta>::iterator;
    using const_iterator = std::vector<SnpMeta>::const_iterator;

    explicit SnpEffects(size_t initial_capacity = 0);

    auto emplace_meta(SnpMeta meta) -> void;

    auto emplace_effects(double additive_effect, double A1_frequency) -> void
    {
        assert(dominance_data_.empty() && "Error: Mixed usage! Dominance mode is active but 2-arg method called.");

        assert(
            additive_data_.size() < snp_meta_.size()
            && "Error: Emplacing more effects than metadata records!");

        additive_data_.push_back(additive_effect);
        frequencies_data_.push_back(A1_frequency);
    }

    auto emplace_effects(
        double additive_effect,
        double dominance_effect,
        double A1_frequency) -> void
    {
        assert((additive_data_.empty() || !dominance_data_.empty()) &&
               "Error: Mixed usage! Additive-only mode is active but 3-arg method called.");

        assert(
            additive_data_.size() < snp_meta_.size()
            && "Error: Emplacing more effects than metadata records!");

        additive_data_.push_back(additive_effect);
        dominance_data_.push_back(dominance_effect);
        frequencies_data_.push_back(A1_frequency);
    }

    auto additive_effects() -> VectorMap
    {
        return VectorMap(
            additive_data_.data(),
            static_cast<Eigen::Index>(additive_data_.size()));
    }
    auto additive_effects() const -> ConstVectorMap
    {
        return ConstVectorMap(
            additive_data_.data(),
            static_cast<Eigen::Index>(additive_data_.size()));
    }

    auto dominance_effects() -> VectorMap
    {
        return VectorMap(
            dominance_data_.data(),
            static_cast<Eigen::Index>(dominance_data_.size()));
    }
    auto dominance_effects() const -> ConstVectorMap
    {
        return ConstVectorMap(
            dominance_data_.data(),
            static_cast<Eigen::Index>(dominance_data_.size()));
    }

    auto frequencies() -> VectorMap
    {
        return VectorMap(
            frequencies_data_.data(),
            static_cast<Eigen::Index>(frequencies_data_.size()));
    }
    auto frequencies() const -> ConstVectorMap
    {
        return ConstVectorMap(
            frequencies_data_.data(),
            static_cast<Eigen::Index>(frequencies_data_.size()));
    }

    auto operator[](size_t index) -> SnpMeta& { return snp_meta_[index]; }
    auto operator[](size_t index) const -> const SnpMeta&
    {
        return snp_meta_[index];
    }

    auto operator[](std::string_view snp_id) -> SnpMeta*;
    auto operator[](std::string_view snp_id) const -> const SnpMeta*;

    auto find_index(std::string_view snp_id) const
        -> std::optional<Eigen::Index>;

    auto shrink_to_fit() -> void;
    auto clear() -> void;

    [[nodiscard]] auto size() const -> size_t { return snp_meta_.size(); }

    auto begin() { return snp_meta_.begin(); }
    auto end() { return snp_meta_.end(); }
    auto begin() const { return snp_meta_.begin(); }
    auto end() const { return snp_meta_.end(); }

   private:
    std::vector<SnpMeta> snp_meta_;
    std::unordered_map<std::string, Eigen::Index> snp_index_map_;

    std::vector<double> additive_data_;
    std::vector<double> dominance_data_;
    std::vector<double> frequencies_data_;
};

}  // namespace gelex
#endif  // GELEX_TYPES_SNP_INFO_H
