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

#ifndef GELEX_DATA_BED_PIPE_VARIANT_DECODER_H_
#define GELEX_DATA_BED_PIPE_VARIANT_DECODER_H_

#include <span>
#include <vector>

#include <Eigen/Core>

namespace gelex::detail
{

class BedVariantDecoder
{
   public:
    BedVariantDecoder(
        Eigen::Index num_raw_samples,
        Eigen::Index bytes_per_variant,
        const std::vector<Eigen::Index>& raw_to_target_sample_idx,
        bool is_dense_mapping);

    void decode(const uint8_t* data_ptr, std::span<double> target_buf) const;

   private:
    void decode_dense(const uint8_t* data_ptr, std::span<double> target_buf)
        const;

    void decode_sparse(const uint8_t* data_ptr, std::span<double> target_buf)
        const;

    Eigen::Index num_raw_samples_ = 0;
    Eigen::Index bytes_per_variant_ = 0;
    const std::vector<Eigen::Index>& raw_to_target_sample_idx_;
    bool is_dense_mapping_ = false;
};

}  // namespace gelex::detail

#endif  // GELEX_DATA_BED_PIPE_VARIANT_DECODER_H_
