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

#ifndef GELEX_DATA_DETAIL_INDEX_PROJECTION_H_
#define GELEX_DATA_DETAIL_INDEX_PROJECTION_H_

#include <Eigen/Core>
#include <span>
#include <string>
#include <vector>

#include "gelex/data/dataframe/index.h"

namespace gelex::detail
{
class IndexProjection
{
   public:
    using index_type = Eigen::Index;
    static constexpr index_type npos = -1;

    IndexProjection(
        const DataFrameIndex<std::string>& source_index,
        const DataFrameIndex<std::string>& target_index);

    IndexProjection(const IndexProjection&) = delete;
    IndexProjection& operator=(const IndexProjection&) = delete;
    IndexProjection(IndexProjection&&) noexcept = default;
    IndexProjection& operator=(IndexProjection&&) noexcept = default;
    ~IndexProjection() = default;

    [[nodiscard]] auto source_size() const -> index_type
    {
        return source_size_;
    }

    [[nodiscard]] auto target_size() const -> index_type
    {
        return static_cast<index_type>(target_to_source_.size());
    }

    [[nodiscard]] auto target_to_source() const -> std::span<const index_type>
    {
        return target_to_source_;
    }

   private:
    index_type source_size_ = 0;
    std::vector<index_type> target_to_source_;
};

}  // namespace gelex::detail

#endif  // GELEX_DATA_DETAIL_INDEX_PROJECTION_H_
