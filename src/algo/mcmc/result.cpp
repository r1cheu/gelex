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

#include "gelex/algo/mcmc/result.h"

#include <cstddef>
#include <ranges>
#include <string>
#include <string_view>
#include <utility>
#include <variant>
#include <vector>

#include <Eigen/Core>

#include "algo/mcmc/detail/pip_computer.h"
#include "algo/mcmc/detail/pve_computer.h"
#include "gelex/algo/mcmc/records.h"
#include "gelex/bayes/model.h"
#include "gelex/exception.h"
#include "gelex/types/genetic_mode.h"

namespace gelex
{

Result::Result(
    Records&& records,
    const BayesModel& model,
    Eigen::Index samples_collected)
    : phenotype_var_(model.phenotype_variance()),
      samples_collected_(samples_collected)
{
    records_ = std::move(records).take_results();
    const auto n_records = records_.size();
    append_derived_records(model, n_records);
    index_records();
}

auto Result::append_derived_records(
    const BayesModel& model,
    std::size_t n_records) -> void
{
    append_pve_records(model, n_records);
    append_pip_records(n_records);
}

auto Result::append_pve_records(const BayesModel& model, std::size_t n_records)
    -> void
{
    if (phenotype_var_ <= 0.0)
    {
        return;
    }

    Eigen::VectorXd additive_beta;
    Eigen::VectorXd dominance_beta;
    const auto pve_computer = detail::PveComputer{model, phenotype_var_};
    for (std::size_t i = 0; i < n_records; ++i)
    {
        const auto& record = records_[i];
        const std::string_view path{record.path};
        if (!std::holds_alternative<RunningStatsResult>(record.value)
            || !path.ends_with("/genetic/coeffs"))
        {
            continue;
        }

        bool is_additive{};
        if (path.contains("/A/"))
        {
            is_additive = true;
        }
        else if (!path.contains("/D/"))
        {
            continue;
        }

        const auto mode = is_additive ? GeneticMode::A : GeneticMode::D;
        const auto& coeffs = std::get<RunningStatsResult>(record.value);
        auto pve = pve_computer.single(mode, coeffs.mean);
        if (mode == GeneticMode::A)
        {
            additive_beta = coeffs.mean;
        }
        else
        {
            dominance_beta = coeffs.mean;
        }

        auto pve_path = std::string{path.substr(0, path.size() - 7)};
        pve_path += "/pve";
        append_record(std::move(pve_path), std::move(pve));
    }
    if (additive_beta.size() != 0 && dominance_beta.size() != 0)
    {
        auto total_pve = pve_computer.total(additive_beta, dominance_beta);
        append_record("state/genetic/pve", std::move(total_pve));
    }
}

auto Result::append_pip_records(std::size_t n_records) -> void
{
    for (std::size_t i = 0; i < n_records; ++i)
    {
        make_pip_records(records_[i]);
    }
}

auto Result::append_record(std::string path, Eigen::VectorXd&& value) -> void
{
    auto size = value.size();
    records_.push_back(
        RecordEntry{
            std::move(path),
            RunningStatsResult{std::move(value), Eigen::VectorXd::Zero(size)},
            std::nullopt});
}

auto Result::append_single_pip_record(
    std::string path,
    const CategoryProbResult& probabilities) -> void
{
    const auto pip_computer = detail::PipComputer{};
    auto pip = pip_computer.single(probabilities.value);
    path += "/pip";
    append_record(std::move(path), std::move(pip));
}

auto Result::append_joint_pip_records(
    std::string path,
    const CategoryProbResult& probabilities) -> void
{
    if (probabilities.value.cols() < 4)
    {
        return;
    }

    const auto pip_computer = detail::PipComputer{};
    auto [additive_pip, dominance_pip]
        = pip_computer.joint(probabilities.value);
    auto additive_pip_path = path;
    additive_pip_path += "/A/pip";
    path += "/D/pip";

    append_record(std::move(additive_pip_path), std::move(additive_pip));
    append_record(std::move(path), std::move(dominance_pip));
}

auto Result::index_records() -> void
{
    record_indices_.reserve(records_.size());
    for (auto [i, record] : std::views::enumerate(records_))
    {
        record_indices_.emplace(record.path, static_cast<std::size_t>(i));
    }
}

auto Result::make_pip_records(const RecordEntry& record) -> void
{
    if (!std::holds_alternative<CategoryProbResult>(record.value))
    {
        return;
    }

    constexpr std::string_view assignment_suffix{"/assignment"};
    const std::string_view path{record.path};
    if (!path.ends_with(assignment_suffix))
    {
        return;
    }

    const auto& probabilities = std::get<CategoryProbResult>(record.value);
    if (probabilities.value.cols() == 0)
    {
        return;
    }

    auto pip_path
        = std::string{path.substr(0, path.size() - assignment_suffix.size())};
    if (path.contains("/joint/"))
    {
        append_joint_pip_records(std::move(pip_path), probabilities);
        return;
    }

    append_single_pip_record(std::move(pip_path), probabilities);
}

auto Result::get(std::string_view path) const -> const RecordResult&
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
