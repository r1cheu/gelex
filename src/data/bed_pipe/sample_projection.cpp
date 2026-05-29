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

#include "sample_projection.h"
#include <Eigen/Core>
#include <cstddef>
#include <string>
#include <vector>
#include "gelex/data/dataframe/index.h"

namespace gelex::genotype::detail
{

SampleProjection::SampleProjection(
    const std::vector<std::string>& raw_ids,
    const dataframe::Index<std::string>& target_index)
{
    const auto num_raw_samples = static_cast<Eigen::Index>(raw_ids.size());

    raw_to_target_sample_idx_.assign(num_raw_samples, -1);

    Eigen::Index mapped_count = 0;
    bool sequential = true;
    for (Eigen::Index i = 0; i < num_raw_samples; ++i)
    {
        const auto& id = raw_ids[static_cast<size_t>(i)];
        if (target_index.contains(id))
        {
            auto target_pos = static_cast<Eigen::Index>(target_index.at(id));
            raw_to_target_sample_idx_[i] = target_pos;
            ++mapped_count;

            if (target_pos != i)
            {
                sequential = false;
            }
        }
        else
        {
            sequential = false;
        }
    }

    num_output_samples_ = static_cast<Eigen::Index>(target_index.size());
    detect_dense_mapping(
        mapped_count,
        sequential
            && (static_cast<size_t>(num_raw_samples) == target_index.size()));
}

auto SampleProjection::detect_dense_mapping(
    Eigen::Index mapped_count,
    bool sequential) -> void
{
    const auto num_raw_samples
        = static_cast<Eigen::Index>(raw_to_target_sample_idx_.size());
    is_dense_mapping_ = sequential && (mapped_count == num_raw_samples);
}

auto SampleProjection::mapping() const -> const std::vector<Eigen::Index>&
{
    return raw_to_target_sample_idx_;
}

auto SampleProjection::is_dense() const -> bool
{
    return is_dense_mapping_;
}

auto SampleProjection::num_output_samples() const -> Eigen::Index
{
    return num_output_samples_;
}

}  // namespace gelex::genotype::detail
