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

#ifndef GELEX_DATA_LOADER_PHENOTYPE_LOADER_H
#define GELEX_DATA_LOADER_PHENOTYPE_LOADER_H

#include <filesystem>
#include <fstream>
#include <string>
#include <unordered_map>
#include <vector>

#include <Eigen/Dense>

namespace gelex::detail
{

class PhenotypeLoader
{
   public:
    PhenotypeLoader(
        const std::filesystem::path& path,
        int pheno_column,
        bool iid_only);

    [[nodiscard]] auto load(
        const std::unordered_map<std::string, Eigen::Index>& id_map) const
        -> Eigen::VectorXd;

    auto name() const -> const std::string& { return name_; }

    auto sample_ids() const -> const std::vector<std::string>&
    {
        return sample_ids_;
    }

   private:
    std::string name_;
    std::vector<std::string> sample_ids_;
    std::vector<double> values_;

    auto init_column(std::ifstream& file, int pheno_column) -> void;
    auto fill_data(std::ifstream& file, int pheno_column, bool iid_only)
        -> void;
};

}  // namespace gelex::detail

#endif  // GELEX_DATA_LOADER_PHENOTYPE_LOADER_H
