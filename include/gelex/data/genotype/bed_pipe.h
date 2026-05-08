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

#ifndef GELEX_DATA_BED_PIPE_H
#define GELEX_DATA_BED_PIPE_H

#include <memory>
#include <string>

#include <Eigen/Core>

#include "gelex/data/dataframe/index.h"

namespace gelex::genotype::detail
{

class SampleProjection;
class BedMmapReader;
class BedVariantDecoder;

}  // namespace gelex::genotype::detail

namespace gelex::genotype
{

class BedPipe
{
   public:
    BedPipe(
        const std::string& bfile_prefix,
        const dataframe::Index<std::string>& sample_index);

    BedPipe(const BedPipe&) = delete;
    BedPipe& operator=(const BedPipe&) = delete;
    BedPipe(BedPipe&&) noexcept;
    BedPipe& operator=(BedPipe&&) noexcept;
    ~BedPipe();

    [[nodiscard]] Eigen::MatrixXd load() const;

    [[nodiscard]] Eigen::MatrixXd load_chunk(
        Eigen::Index start_col,
        Eigen::Index end_col) const;

    void load_chunk(
        Eigen::Ref<Eigen::MatrixXd> target_buf,
        Eigen::Index start_col,
        Eigen::Index end_col) const;

    [[nodiscard]] Eigen::MatrixXd select(
        std::span<const Eigen::Index> col_indices) const;

    [[nodiscard]] Eigen::Index num_samples() const;
    [[nodiscard]] Eigen::Index num_snps() const;

   private:
    std::unique_ptr<gelex::genotype::detail::SampleProjection> projection_;
    std::unique_ptr<gelex::genotype::detail::BedMmapReader> bed_reader_;
    std::unique_ptr<gelex::genotype::detail::BedVariantDecoder> decoder_;
    Eigen::Index num_raw_snps_ = 0;
    Eigen::Index num_output_samples_ = 0;
};

}  // namespace gelex::genotype
#endif  // GELEX_DATA_BED_PIPE_H
