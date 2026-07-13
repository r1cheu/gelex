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

#ifndef APPS_CLI_REML_DATA_H_
#define APPS_CLI_REML_DATA_H_

#include <Eigen/Core>
#include <map>
#include <optional>
#include <string>
#include <vector>

#include "gelex/data/dataframe/dataframe.h"
#include "gelex/data/dataframe/index.h"
#include "gelex/freq/design.h"

namespace cli
{

struct RemlDataConfig
{
    std::vector<std::string> grm;
    std::optional<std::string> drand_path;
    std::vector<std::string> qrand_paths;
    std::vector<std::string> interactions;  // "<name_a>:<name_b>" pairs
};

// Loads the random-effect data a REML fit consumes (discrete/quantitative/GRM)
// into design matrices. Satisfies the BaseDataHandler concept so it can drive
// load_base_data directly; compose it when a command needs extra data sources.
class RemlDataLoader
{
   public:
    explicit RemlDataLoader(const RemlDataConfig& config) noexcept;

    auto load_indices(
        std::vector<const gelex::DataFrameIndex<std::string>*>& indices)
        -> void;
    auto gather(const gelex::DataFrameIndex<std::string>& common_index) -> void;
    auto results() && -> std::vector<gelex::freq::RandomDesign>;

    auto grm_indices() const noexcept
        -> const std::vector<gelex::DataFrameIndex<std::string>>&
    {
        return grm_indices_;
    }

   private:
    // GRM referenced only by an interaction operand: it joins the sample-index
    // intersection like --grm but is consumed by the interaction rather than
    // emitted as a standalone variance component. Keyed by GRM prefix.
    struct InteractionGrm
    {
        gelex::DataFrameIndex<std::string> index;
        Eigen::MatrixXd K;
    };

    // Resolves an interaction operand to its kernel: a path operand from
    // interaction_grms_, otherwise a loaded effect matched by name in
    // random_designs_. Throws listing loaded effects when neither matches.
    auto resolve_operand(const std::string& name) const
        -> const Eigen::MatrixXd&;

    const RemlDataConfig& config_;
    std::vector<gelex::DataFrameIndex<std::string>> grm_indices_;
    std::vector<gelex::freq::RandomDesign> random_designs_;
    std::optional<gelex::DataFrame<std::string>> drand_;
    std::vector<gelex::DataFrame<std::string>> qrand_;
    std::map<std::string, InteractionGrm> interaction_grms_;
};

}  // namespace cli

#endif  // APPS_CLI_REML_DATA_H_
