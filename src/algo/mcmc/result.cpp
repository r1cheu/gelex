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
#include <optional>
#include <ranges>
#include <string>
#include <string_view>
#include <utility>
#include <variant>
#include <vector>

#include <Eigen/Core>

#include "gelex/algo/mcmc/records.h"
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
    append_derived_records(model, n_records);
    index_records();
}

auto mcmc::Result::append_derived_records(
    const BayesModel& model,
    std::size_t n_records) -> void
{
    bool has_total{};
    for (std::size_t i = 0; i < n_records; ++i)
    {
        auto pve_data = make_pve(model, records_[i]);
        if (pve_data && phenotype_var_ > 0.0)
        {
            auto [genetic_values, pve_path] = std::move(*pve_data);
            append_pve_record(std::move(pve_path), genetic_values);

            if (pve_buffer_.size() == 0)
            {
                pve_buffer_ = std::move(genetic_values);
            }
            else
            {
                has_total = true;
                pve_buffer_ += genetic_values;
            }
        }
        make_pip_records(records_[i]);
    }
    if (has_total)
    {
        append_pve_record("state/genetic/pve", pve_buffer_);
    }
}

auto mcmc::Result::append_record(std::string path, Eigen::VectorXd&& value)
    -> void
{
    auto size = value.size();
    records_.push_back(
        RecordEntry{
            std::move(path),
            stats::RunningStatsResult{
                std::move(value), Eigen::VectorXd::Zero(size)},
            std::nullopt});
}

auto mcmc::Result::append_pve_record(
    std::string path,
    const Eigen::Ref<const Eigen::MatrixXd>& genetic_values) -> void
{
    auto pve = (stats::detail::matvar(
                    genetic_values, stats::detail::VarNormType::Population)
                    .transpose() /= phenotype_var_)
                   .eval();
    append_record(std::move(path), std::move(pve));
}

auto mcmc::Result::index_records() -> void
{
    record_indices_.reserve(records_.size());
    for (auto [i, record] : std::views::enumerate(records_))
    {
        record_indices_.emplace(record.path, static_cast<std::size_t>(i));
    }
}

auto mcmc::Result::make_pve(const BayesModel& model, const RecordEntry& record)
    -> std::optional<std::pair<Eigen::MatrixXd, std::string>>
{
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
    if (design->X.cols() != coeffs.mean.size())
    {
        throw GelexException(
            "Result: genetic coefficient size does not match design");
    }

    // Drop the trailing "/coeffs" segment.
    auto pve_path = std::string{path.substr(0, path.size() - 7)};
    pve_path += "/pve";
    return std::pair{
        (design->X.matrix() * coeffs.mean.asDiagonal()).eval(),
        std::move(pve_path)};
}

auto mcmc::Result::make_pip_records(const RecordEntry& record) -> void
{
    if (!std::holds_alternative<stats::CategoryProbResult>(record.value))
    {
        return;
    }

    constexpr std::string_view assignment_suffix{"/assignment"};
    const std::string_view path{record.path};
    if (!path.ends_with(assignment_suffix))
    {
        return;
    }

    const auto& probabilities
        = std::get<stats::CategoryProbResult>(record.value);
    if (probabilities.value.cols() == 0)
    {
        return;
    }

    auto pip_path
        = std::string{path.substr(0, path.size() - assignment_suffix.size())};
    if (path.contains("/joint/"))
    {
        if (probabilities.value.cols() < 4)
        {
            return;
        }
        auto additive_pip
            = (probabilities.value.col(1) + probabilities.value.col(3)).eval();
        auto dominance_pip
            = (probabilities.value.col(2) + probabilities.value.col(3)).eval();
        auto additive_pip_path = pip_path;
        additive_pip_path += "/A/pip";
        pip_path += "/D/pip";

        append_record(std::move(additive_pip_path), std::move(additive_pip));
        append_record(std::move(pip_path), std::move(dominance_pip));
        return;
    }

    auto pip = (Eigen::VectorXd::Ones(probabilities.value.rows())
                - probabilities.value.col(0))
                   .eval();
    pip_path += "/pip";
    append_record(std::move(pip_path), std::move(pip));
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
