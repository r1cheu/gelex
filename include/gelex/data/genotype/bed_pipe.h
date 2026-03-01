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

#include <filesystem>
#include <memory>

#include <Eigen/Core>

#include "gelex/data/genotype/sample_manager.h"

namespace gelex
{

namespace detail
{

class SampleProjection;
class BedMmapReader;
class BedVariantDecoder;

}  // namespace detail

class BedPipe
{
   public:
    BedPipe(
        const std::filesystem::path& bed_prefix,
        std::shared_ptr<SampleManager> sample_manager);

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

    [[nodiscard]] Eigen::Index num_samples() const;
    [[nodiscard]] Eigen::Index num_snps() const;

   private:
    std::shared_ptr<SampleManager> sample_manager_;
    std::unique_ptr<detail::SampleProjection> projection_;
    std::unique_ptr<detail::BedMmapReader> bed_reader_;
    std::unique_ptr<detail::BedVariantDecoder> decoder_;
    Eigen::Index num_raw_snps_ = 0;
};

}  // namespace gelex
#endif  // GELEX_DATA_BED_PIPE_H
