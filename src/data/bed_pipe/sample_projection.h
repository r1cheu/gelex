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

#ifndef GELEX_DATA_BED_PIPE_SAMPLE_PROJECTION_H_
#define GELEX_DATA_BED_PIPE_SAMPLE_PROJECTION_H_

#include <memory>
#include <string>
#include <vector>

#include <Eigen/Core>

#include "gelex/data/sample_manager.h"

namespace gelex::detail
{

class SampleProjection
{
   public:
    SampleProjection(
        const std::vector<std::string>& raw_ids,
        const std::shared_ptr<SampleManager>& sample_manager);

    [[nodiscard]] auto mapping() const -> const std::vector<Eigen::Index>&;

    [[nodiscard]] auto is_dense() const -> bool;

   private:
    std::vector<Eigen::Index> raw_to_target_sample_idx_;
    bool is_dense_mapping_ = false;
};

}  // namespace gelex::detail

#endif  // GELEX_DATA_BED_PIPE_SAMPLE_PROJECTION_H_
