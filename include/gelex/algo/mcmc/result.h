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

#ifndef GELEX_ALGO_MCMC_RESULT_H_
#define GELEX_ALGO_MCMC_RESULT_H_

#include <Eigen/Core>
#include <cstddef>
#include <span>
#include <string>
#include <string_view>
#include <unordered_map>
#include <vector>

#include "gelex/algo/mcmc/records.h"

namespace gelex
{

class BayesModel;

class Result
{
   public:
    Result(
        Records&& records,
        const BayesModel& model,
        Eigen::Index samples_collected);

    auto get(std::string_view path) const -> const RecordResult&;
    auto records() const -> std::span<const RecordEntry> { return records_; }
    auto samples_collected() const -> Eigen::Index
    {
        return samples_collected_;
    }
    auto phenotype_variance() const -> double { return phenotype_var_; }

   private:
    auto append_derived_records(const BayesModel& model, std::size_t n_records)
        -> void;
    auto append_pve_records(const BayesModel& model, std::size_t n_records)
        -> void;
    auto append_pip_records(std::size_t n_records) -> void;
    auto append_record(std::string path, Eigen::VectorXd&& value) -> void;
    auto append_single_pip_record(
        std::string path,
        const CategoryProbResult& probabilities) -> void;
    auto append_joint_pip_records(
        std::string path,
        const CategoryProbResult& probabilities) -> void;
    auto index_records() -> void;

    auto make_pip_records(const RecordEntry& record) -> void;

    double phenotype_var_;

    std::vector<RecordEntry> records_;
    std::unordered_map<std::string, std::size_t> record_indices_;
    Eigen::Index samples_collected_{};
};

}  // namespace gelex

#endif  // GELEX_ALGO_MCMC_RESULT_H_
