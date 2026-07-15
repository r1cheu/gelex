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

#ifndef GELEX_DATA_BED_H_
#define GELEX_DATA_BED_H_

#include <Eigen/Core>
#include <string>
#include <utility>

#include "gelex/data/dataframe/dataframe.h"
#include "gelex/data/dataframe/index.h"
#include "gelex/data/detail/bed_source.h"
#include "gelex/data/detail/index_projection.h"

namespace gelex
{

class LocusEncoder;

// Owns a PLINK1 .bed dataset's 2-bit-packed source, its .fam/.bim metadata, and
// the sample projection onto the current analysis order. Bed exposes only
// metadata and the sample view; genotype access runs through LocusEncoder.
class Bed
{
   public:
    Bed(const Bed&) = delete;
    Bed& operator=(const Bed&) = delete;
    Bed(Bed&&) noexcept = default;
    Bed& operator=(Bed&&) noexcept = default;
    ~Bed() = default;

    [[nodiscard]] auto num_samples() const noexcept -> Eigen::Index
    {
        return index_projection_.target_size();
    }

    [[nodiscard]] auto num_snps() const noexcept -> Eigen::Index
    {
        return static_cast<Eigen::Index>(bed_source_.size());
    }

    [[nodiscard]] auto sample_index() const noexcept
        -> const DataFrameIndex<std::string>&
    {
        return sample_index_;
    }

    [[nodiscard]] auto snp_index() const noexcept
        -> const DataFrameIndex<std::string>&
    {
        return bim_.index();
    }

    [[nodiscard]] auto bim() const noexcept -> const DataFrame<std::string>&
    {
        return bim_;
    }

    // Narrows the sample view to target, rebuilding the source->target
    // projection. target must be a subset of the source .fam samples.
    auto gather(const DataFrameIndex<std::string>& target) -> void
    {
        index_projection_ = detail::IndexProjection{source_index_, target};
        sample_index_ = target;
    }

   private:
    Bed(detail::BedSource bed_source,
        DataFrameIndex<std::string> source_index,
        DataFrame<std::string> bim)
        : bed_source_{std::move(bed_source)},
          source_index_{std::move(source_index)},
          index_projection_{source_index_, source_index_},
          sample_index_{source_index_},
          bim_{std::move(bim)}
    {
    }

    detail::BedSource bed_source_;
    DataFrameIndex<std::string> source_index_;
    detail::IndexProjection index_projection_;
    DataFrameIndex<std::string> sample_index_;
    DataFrame<std::string> bim_;

    friend auto open_bed(const std::string& bfile_prefix) -> Bed;
    friend class LocusEncoder;
};

[[nodiscard]] auto open_bed(const std::string& bfile_prefix) -> Bed;

}  // namespace gelex

#endif  // GELEX_DATA_BED_H_
