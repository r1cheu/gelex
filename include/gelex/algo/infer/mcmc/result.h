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

#ifndef GELEX_ALGO_INFER_MCMC_RESULT_H_
#define GELEX_ALGO_INFER_MCMC_RESULT_H_

#include <cstddef>
#include <string>
#include <string_view>
#include <unordered_map>
#include <vector>

#include <Eigen/Core>

#include "gelex/algo/infer/mcmc/records.h"
#include "gelex/types/genetic_effect_type.h"

namespace gelex
{

class BayesModel;
struct FixedSummary;
struct GeneticSummary;
struct PosteriorSummary;
struct RandomSummary;

namespace mcmc
{

class Result
{
   public:
    Result(
        mcmc::Records&& records,
        const BayesModel& model,
        Eigen::Index samples_collected);

    auto get(std::string_view path) const -> const RecordResult&;
    auto samples_collected() const -> Eigen::Index
    {
        return samples_collected_;
    }
    auto phenotype_variance() const -> double { return phenotype_var_; }
    auto allele_frequency() const -> const Eigen::VectorXd& { return p_freq_; }

    auto fixed() const -> const FixedSummary&;
    auto random() const -> const std::vector<RandomSummary>&;
    auto genetics() const -> const std::vector<GeneticSummary>&;
    auto genetic(GeneticMode type) const -> const GeneticSummary*;
    auto residual() const -> const PosteriorSummary&;
    auto allele_freq() const -> const Eigen::VectorXd&;

   private:
    double phenotype_var_;

    Eigen::VectorXd p_freq_;
    std::vector<RecordEntry> records_;
    std::unordered_map<std::string, std::size_t> record_indices_;
    Eigen::Index samples_collected_{};
};

}  // namespace mcmc

}  // namespace gelex

#endif  // GELEX_ALGO_INFER_MCMC_RESULT_H_
