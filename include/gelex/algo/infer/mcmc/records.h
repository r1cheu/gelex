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

#ifndef GELEX_ALGO_INFER_MCMC_RECORDS_H_
#define GELEX_ALGO_INFER_MCMC_RECORDS_H_

#include <cstddef>
#include <optional>
#include <span>
#include <string>
#include <string_view>
#include <unordered_map>
#include <variant>
#include <vector>

#include <Eigen/Core>

#include "gelex/algo/infer/mcmc/detail/draw_stats.h"
#include "gelex/infra/field_visitor.h"
#include "gelex/infra/stats/result.h"

namespace gelex
{

class BayesState;
class BayesModel;

namespace bayes
{
class BayesPrior;
}

namespace mcmc
{

using RecordResult
    = std::variant<stats::RunningStatsResult, stats::CategoryProbResult>;

struct RecordEntry
{
    std::string path;
    RecordResult value;
    std::optional<std::vector<std::string>> names;
};

class Records : private infra::FieldVisitor
{
   public:
    Records() = default;
    Records(const Records&) = delete;
    auto operator=(const Records&) -> Records& = delete;
    Records(Records&& other) noexcept;
    auto operator=(Records&& other) noexcept -> Records&;
    ~Records() override = default;

    auto store(
        const BayesModel& model,
        const bayes::BayesPrior& prior,
        BayesState& state) -> void;

    auto take_results() && -> std::vector<RecordEntry>;

   private:
    struct ScalarRecord
    {
        detail::ScalarDrawStats draws;
    };

    struct VectorRecord
    {
        detail::VectorDrawStats draws;
    };

    struct CategoricalRecord
    {
        CategoricalRecord(Eigen::Index n_items, Eigen::Index n_categories);

        detail::CategoricalDrawStats draws;
    };

    using Record = std::variant<ScalarRecord, VectorRecord, CategoricalRecord>;

    auto result(std::string_view path) -> RecordResult;

    template <typename RecordType, typename Value>
    auto store_record(std::string_view name, Value&& value) -> void;

    auto on(
        std::string_view name,
        Eigen::Ref<Eigen::VectorXf> value,
        FieldFlag flags) -> void override;
    auto on(
        std::string_view name,
        Eigen::Ref<Eigen::VectorXd> value,
        FieldFlag flags) -> void override;
    auto on(
        std::string_view name,
        Eigen::Ref<Eigen::VectorXi> value,
        FieldFlag flags) -> void override;
    auto on(std::string_view name, double& value, FieldFlag flags)
        -> void override;
    auto on(std::string_view name, int& value, FieldFlag flags)
        -> void override;
    auto on(
        std::string_view name,
        std::span<const std::string> value,
        FieldFlag flags) -> void override;

    std::vector<Record> records_;
    std::vector<std::string> paths_;
    std::vector<std::optional<std::vector<std::string>>> names_;
    std::unordered_map<std::string, std::size_t> indices_;
    std::unordered_map<std::string, Eigen::Index> category_counts_;
};

}  // namespace mcmc

}  // namespace gelex

#endif  // GELEX_ALGO_INFER_MCMC_RECORDS_H_
