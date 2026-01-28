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

#ifndef GELEX_DATA_LOADER_QCOVARIATE_LOADER_H
#define GELEX_DATA_LOADER_QCOVARIATE_LOADER_H

#include <filesystem>
#include <fstream>
#include <string>
#include <unordered_map>
#include <vector>

#include <Eigen/Core>

#include "../src/types/covariates.h"

namespace gelex::detail
{

class QuantitativeCovariateLoader
{
   public:
    QuantitativeCovariateLoader(
        const std::filesystem::path& path,
        bool iid_only);

    [[nodiscard]] auto load(
        const std::unordered_map<std::string, Eigen::Index>& id_map) const
        -> QuantitativeCovariate;

    auto sample_ids() const -> const std::vector<std::string>&
    {
        return sample_ids_;
    }

    auto column_names() const -> const std::vector<std::string>&
    {
        return column_names_;
    }

   private:
    static constexpr size_t IdColumnCount = 2;

    std::vector<std::string> column_names_;
    std::vector<std::string> sample_ids_;
    std::vector<std::vector<double>> columns_;

    auto init_columns(std::ifstream& file) -> void;
    auto fill_columns(std::ifstream& file, bool iid_only) -> void;
    [[nodiscard]] static auto is_valid_covariate_value(double value) -> bool;
};

}  // namespace gelex::detail

#endif  // GELEX_DATA_LOADER_QCOVARIATE_LOADER_H
