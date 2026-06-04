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

#include "gelex/algo/infer/mcmc/result.h"

#include <cstddef>
#include <optional>
#include <ranges>
#include <string>
#include <string_view>
#include <utility>
#include <variant>
#include <vector>

#include <Eigen/Core>

#include "gelex/algo/infer/mcmc/records.h"
#include "gelex/bayes/model.h"
#include "gelex/exception.h"
#include "gelex/infra/stats/detail/var.h"
#include "gelex/types/genetic_effect_type.h"

namespace gelex
{

mcmc::Result::Result(
    mcmc::Records&& records,
    const BayesModel& model,
    Eigen::Index samples_collected)
    : phenotype_var_(model.phenotype_variance()),
      samples_collected_(samples_collected)
{
    records_ = std::move(records).take_results();
    const auto n_records = records_.size();
    records_.reserve(n_records * 3);
    for (std::size_t i = 0; i < n_records; ++i)
    {
        auto pve_record = make_pve_record(model, records_[i], phenotype_var_);
        if (pve_record)
        {
            records_.push_back(std::move(*pve_record));
        }
        auto pip_records = make_pip_records(records_[i]);
        for (auto& pip_record : pip_records)
        {
            records_.push_back(std::move(pip_record));
        }
    }

    record_indices_.reserve(records_.size());
    for (auto [i, record] : std::views::enumerate(records_))
    {
        record_indices_.emplace(record.path, static_cast<std::size_t>(i));
    }
}

auto mcmc::Result::make_pve_record(
    const BayesModel& model,
    const RecordEntry& record,
    double phenotype_var) -> std::optional<RecordEntry>
{
    if (phenotype_var <= 0.0)
    {
        return std::nullopt;
    }
    if (!std::holds_alternative<stats::RunningStatsResult>(record.value))
    {
        return std::nullopt;
    }

    const std::string_view path{record.path};
    if (!path.ends_with("/genetic/coeffs"))
    {
        return std::nullopt;
    }

    const bayes::GeneticDesign* design{};
    if (path.contains("/A/"))
    {
        design = model.genetic(GeneticMode::A);
    }
    else if (path.contains("/D/"))
    {
        design = model.genetic(GeneticMode::D);
    }
    if (design == nullptr)
    {
        return std::nullopt;
    }

    const auto& coeffs = std::get<stats::RunningStatsResult>(record.value);
    auto pve
        = stats::detail::matvar(design->X.matrix() * coeffs.mean.asDiagonal())
              .transpose()
              .eval();
    pve.array() /= phenotype_var;

    // Drop the trailing "/coeffs" segment.
    auto pve_path = std::string{path.substr(0, path.size() - 7)};
    pve_path += "/pve";
    return RecordEntry{
        std::move(pve_path),
        stats::RunningStatsResult{
            std::move(pve), Eigen::VectorXd::Zero(coeffs.mean.size())},
        std::nullopt};
}

auto mcmc::Result::make_pip_records(const RecordEntry& record)
    -> std::vector<RecordEntry>
{
    if (!std::holds_alternative<stats::CategoryProbResult>(record.value))
    {
        return {};
    }

    constexpr std::string_view assignment_suffix{"/assignment"};
    const std::string_view path{record.path};
    if (!path.ends_with(assignment_suffix))
    {
        return {};
    }

    const auto& probabilities
        = std::get<stats::CategoryProbResult>(record.value);
    if (probabilities.value.cols() == 0)
    {
        return {};
    }

    auto pip_path
        = std::string{path.substr(0, path.size() - assignment_suffix.size())};
    std::vector<RecordEntry> output;
    if (path.contains("/joint/"))
    {
        if (probabilities.value.cols() < 4)
        {
            return {};
        }

        output.reserve(2);
        auto additive_pip
            = (probabilities.value.col(1) + probabilities.value.col(3)).eval();
        auto additive_pip_path = pip_path;
        additive_pip_path += "/A/pip";
        output.push_back(
            RecordEntry{
                std::move(additive_pip_path),
                stats::RunningStatsResult{
                    std::move(additive_pip),
                    Eigen::VectorXd::Zero(probabilities.value.rows())},
                std::nullopt});

        auto dominance_pip
            = (probabilities.value.col(2) + probabilities.value.col(3)).eval();
        pip_path += "/D/pip";
        output.push_back(
            RecordEntry{
                std::move(pip_path),
                stats::RunningStatsResult{
                    std::move(dominance_pip),
                    Eigen::VectorXd::Zero(probabilities.value.rows())},
                std::nullopt});
        return output;
    }

    auto pip = (Eigen::VectorXd::Ones(probabilities.value.rows())
                - probabilities.value.col(0))
                   .eval();

    pip_path += "/pip";
    output.push_back(
        RecordEntry{
            std::move(pip_path),
            stats::RunningStatsResult{
                std::move(pip),
                Eigen::VectorXd::Zero(probabilities.value.rows())},
            std::nullopt});
    return output;
}

auto mcmc::Result::get(std::string_view path) const -> const RecordResult&
{
    auto key = std::string{path};
    const auto it = record_indices_.find(key);
    if (it == record_indices_.end())
    {
        throw GelexException("Result: missing record " + key);
    }
    return records_[it->second].value;
}

}  // namespace gelex
