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

#ifndef GELEX_DATA_LOADER_SAMPLE_LOADER_H
#define GELEX_DATA_LOADER_SAMPLE_LOADER_H

#include <filesystem>
#include <string>
#include <unordered_map>
#include <vector>

#include <Eigen/Dense>

namespace gelex::detail
{

// -----------------------------------------------------------------------------
// FamLoader (Sample Info)
// -----------------------------------------------------------------------------
class FamLoader
{
   public:
    explicit FamLoader(const std::filesystem::path& path, bool iid_only);

    const std::vector<std::string>& ids() const { return ids_; }
    const std::unordered_map<std::string, Eigen::Index>& data() const
    {
        return data_;
    }
    std::vector<std::string>&& take_ids() && { return std::move(ids_); }

   private:
    void set_ids(const std::filesystem::path& path, bool iid_only);
    void set_index_map();
    std::vector<std::string> ids_;
    std::unordered_map<std::string, Eigen::Index> data_;
};

}  // namespace gelex::detail

#endif  // GELEX_DATA_LOADER_SAMPLE_LOADER_H
