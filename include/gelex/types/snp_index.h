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

#ifndef GELEX_TYPES_SNP_INDEX_H
#define GELEX_TYPES_SNP_INDEX_H

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

class SnpIndex
{
   public:
    using iterator = std::vector<SnpMeta>::iterator;
    using const_iterator = std::vector<SnpMeta>::const_iterator;

    explicit SnpIndex(size_t initial_capacity = 0);

    auto emplace(SnpMeta meta) -> void;

    auto operator[](size_t index) -> SnpMeta& { return snp_meta_[index]; }
    auto operator[](size_t index) const -> const SnpMeta&
    {
        return snp_meta_[index];
    }

    auto find(std::string_view snp_id) -> SnpMeta*;
    auto find(std::string_view snp_id) const -> const SnpMeta*;

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
};

}  // namespace gelex
#endif  // GELEX_TYPES_SNP_INDEX_H
