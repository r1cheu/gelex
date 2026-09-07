// Copyright 2026 RuLei Chen
// SPDX-License-Identifier: Apache-2.0

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
