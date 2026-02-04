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

#ifndef GELEX_DATA_SAMPLE_MANAGER_H_
#define GELEX_DATA_SAMPLE_MANAGER_H_

#include <filesystem>
#include <memory>
#include <span>
#include <string>
#include <unordered_map>
#include <vector>

#include <Eigen/Core>

namespace gelex
{

class SampleManager
{
   public:
    explicit SampleManager(
        const std::filesystem::path& fam_path,
        bool iid_only = false);

    SampleManager(const SampleManager&) = delete;
    SampleManager(SampleManager&&) noexcept = default;
    SampleManager& operator=(const SampleManager&) = delete;
    SampleManager& operator=(SampleManager&&) noexcept = default;
    ~SampleManager() = default;

    void intersect(std::span<const std::string> ids);

    void finalize();

    static auto create_finalized(
        const std::filesystem::path& bed_path,
        bool iid_only = false) -> std::shared_ptr<SampleManager>;

    [[nodiscard]] const std::vector<std::string>& common_ids() const
    {
        return common_ids_;
    }

    [[nodiscard]] const std::unordered_map<std::string, Eigen::Index>&
    common_id_map() const
    {
        return common_id_map_;
    }

    [[nodiscard]] size_t num_common_samples() const
    {
        return common_ids_.size();
    }
    [[nodiscard]] bool has_common_samples() const
    {
        return !common_ids_.empty();
    }

   private:
    std::vector<std::string> common_ids_;
    std::unordered_map<std::string, Eigen::Index> common_id_map_;
};

}  // namespace gelex

#endif  // GELEX_DATA_SAMPLE_MANAGER_H_
