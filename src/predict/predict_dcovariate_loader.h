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

#ifndef GELEX_PREDICT_PREDICT_DCOVARIATE_LOADER_H
#define GELEX_PREDICT_PREDICT_DCOVARIATE_LOADER_H

#include <filesystem>
#include <map>
#include <string>
#include <unordered_map>
#include <vector>

#include <Eigen/Dense>

namespace gelex::detail
{

class DcovarPredictLoader
{
   public:
    explicit DcovarPredictLoader(
        const std::filesystem::path& path,
        bool iid_only);

    auto load(const std::unordered_map<std::string, Eigen::Index>& id_map) const
        -> std::map<std::string, std::vector<std::string>>;

    const std::vector<std::string>& names() const { return names_; }
    const std::unordered_map<std::string, std::vector<std::string>>& data()
        const
    {
        return data_;
    }

   private:
    void set_names(std::ifstream& file);
    void set_data(std::ifstream& file, bool iid_only);

    std::vector<std::string> names_;
    std::unordered_map<std::string, std::vector<std::string>> data_;
};
}  // namespace gelex::detail

#endif  // GELEX_PREDICT_PREDICT_DCOVARIATE_LOADER_H
