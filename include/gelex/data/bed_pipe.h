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
#include <span>
#include <vector>

#include <Eigen/Core>

#include "gelex/data/sample_manager.h"
#include "mio.h"

namespace gelex
{

class BedPipe
{
   public:
    BedPipe(
        const std::filesystem::path& bed_prefix,
        std::shared_ptr<SampleManager> sample_manager);

    BedPipe(const BedPipe&) = delete;
    BedPipe& operator=(const BedPipe&) = delete;
    BedPipe(BedPipe&&) noexcept = default;
    BedPipe& operator=(BedPipe&&) noexcept = default;
    ~BedPipe() = default;

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

    static auto format_bed_path(std::string_view bed_path)
        -> std::filesystem::path;

   private:
    mio::mmap_source mmap_;
    std::shared_ptr<SampleManager> sample_manager_;

    std::vector<Eigen::Index> raw_to_target_sample_idx_;

    bool is_dense_mapping_ = false;

    Eigen::Index num_raw_samples_ = 0;
    Eigen::Index num_raw_snps_ = 0;
    Eigen::Index bytes_per_variant_ = 0;
    std::filesystem::path bed_path_;

    void decode_variant_dense(
        const uint8_t* data_ptr,
        std::span<double> target_buf) const;

    void decode_variant_sparse(
        const uint8_t* data_ptr,
        std::span<double> target_buf) const;

    void init_bed_mmap(const std::filesystem::path& bed_path);
};

}  // namespace gelex
#endif  // GELEX_DATA_BED_PIPE_H
