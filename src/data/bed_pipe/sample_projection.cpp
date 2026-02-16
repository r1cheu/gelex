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

namespace gelex::detail
{

SampleProjection::SampleProjection(
    const std::vector<std::string>& raw_ids,
    const std::shared_ptr<SampleManager>& sample_manager)
{
    const auto num_raw_samples = static_cast<Eigen::Index>(raw_ids.size());

    raw_to_target_sample_idx_.assign(num_raw_samples, -1);

    const auto& target_map = sample_manager->common_id_map();

    Eigen::Index mapped_count = 0;
    bool sequential = true;
    for (Eigen::Index i = 0; i < num_raw_samples; ++i)
    {
        if (auto it = target_map.find(raw_ids[static_cast<size_t>(i)]);
            it != target_map.end())
        {
            raw_to_target_sample_idx_[i] = it->second;
            ++mapped_count;

            if (it->second != i)
            {
                sequential = false;
            }
        }
        else
        {
            sequential = false;
        }
    }

    is_dense_mapping_ = sequential && (mapped_count == num_raw_samples)
                        && (static_cast<size_t>(num_raw_samples)
                            == sample_manager->num_common_samples());
}

auto SampleProjection::mapping() const -> const std::vector<Eigen::Index>&
{
    return raw_to_target_sample_idx_;
}

auto SampleProjection::is_dense() const -> bool
{
    return is_dense_mapping_;
}

}  // namespace gelex::detail
