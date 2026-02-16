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

#include "variant_decoder.h"

#include <cstring>

#include "data/decode_lut.h"

namespace gelex::detail
{

BedVariantDecoder::BedVariantDecoder(
    Eigen::Index num_raw_samples,
    Eigen::Index bytes_per_variant,
    const std::vector<Eigen::Index>& raw_to_target_sample_idx,
    bool is_dense_mapping)
    : num_raw_samples_(num_raw_samples),
      bytes_per_variant_(bytes_per_variant),
      raw_to_target_sample_idx_(raw_to_target_sample_idx),
      is_dense_mapping_(is_dense_mapping)
{
}

void BedVariantDecoder::decode(
    const uint8_t* data_ptr,
    std::span<double> target_buf) const
{
    if (is_dense_mapping_)
    {
        decode_dense(data_ptr, target_buf);
        return;
    }

    decode_sparse(data_ptr, target_buf);
}

void BedVariantDecoder::decode_dense(
    const uint8_t* data_ptr,
    std::span<double> target_buf) const
{
    const auto& lut = kDecodeLut;
    const Eigen::Index num_bytes = bytes_per_variant_;
    const Eigen::Index total_samples = num_raw_samples_;

    for (Eigen::Index i = 0; i < num_bytes; ++i)
    {
        const uint8_t byte = data_ptr[i];
        const auto& vals = lut[byte];

        const Eigen::Index base_idx = i * 4;
        if (base_idx + 4 <= total_samples)
        {
            std::memcpy(&target_buf[base_idx], vals.data(), 4 * sizeof(double));
        }
        else
        {
            for (int k = 0; k < 4 && (base_idx + k) < total_samples; ++k)
            {
                target_buf[base_idx + k] = vals[k];
            }
        }
    }
}

void BedVariantDecoder::decode_sparse(
    const uint8_t* data_ptr,
    std::span<double> target_buf) const
{
    const auto& lut = kDecodeLut;
    const Eigen::Index num_bytes = bytes_per_variant_;

    for (Eigen::Index i = 0; i < num_bytes; ++i)
    {
        const uint8_t byte = data_ptr[i];
        const auto& vals = lut[byte];

        const Eigen::Index base_raw_idx = i * 4;

        if (Eigen::Index tidx = raw_to_target_sample_idx_[base_raw_idx];
            tidx != -1)
        {
            target_buf[tidx] = vals[0];
        }

        if (base_raw_idx + 1 < num_raw_samples_) [[likely]]
        {
            if (Eigen::Index tidx = raw_to_target_sample_idx_[base_raw_idx + 1];
                tidx != -1)
            {
                target_buf[tidx] = vals[1];
            }
        }

        if (base_raw_idx + 2 < num_raw_samples_) [[likely]]
        {
            if (Eigen::Index tidx = raw_to_target_sample_idx_[base_raw_idx + 2];
                tidx != -1)
            {
                target_buf[tidx] = vals[2];
            }
        }

        if (base_raw_idx + 3 < num_raw_samples_) [[likely]]
        {
            if (Eigen::Index tidx = raw_to_target_sample_idx_[base_raw_idx + 3];
                tidx != -1)
            {
                target_buf[tidx] = vals[3];
            }
        }
    }
}

}  // namespace gelex::detail
