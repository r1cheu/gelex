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
#include <ranges>
#include <string>
#include <string_view>
#include <utility>

#include <Eigen/Core>

#include "gelex/algo/infer/mcmc/records.h"
#include "gelex/bayes/model.h"
#include "gelex/exception.h"
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
    if (const auto* design = model.genetic(GeneticMode::A); design)
    {
        p_freq_ = design->X.mean().array() / 2;
    }

    records_ = records.take_results();
    record_indices_.reserve(records_.size());
    for (auto [i, record] : std::views::enumerate(records_))
    {
        record_indices_.emplace(record.path, static_cast<std::size_t>(i));
    }
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

auto mcmc::Result::fixed() const -> const FixedSummary&
{
    throw GelexException("Result::fixed is legacy and will be removed");
}

auto mcmc::Result::random() const -> const std::vector<RandomSummary>&
{
    throw GelexException("Result::random is legacy and will be removed");
}

auto mcmc::Result::genetics() const -> const std::vector<GeneticSummary>&
{
    throw GelexException("Result::genetics is legacy and will be removed");
}

auto mcmc::Result::genetic(GeneticMode type) const -> const GeneticSummary*
{
    static_cast<void>(type);
    throw GelexException("Result::genetic is legacy and will be removed");
}

auto mcmc::Result::residual() const -> const PosteriorSummary&
{
    throw GelexException("Result::residual is legacy and will be removed");
}

auto mcmc::Result::allele_freq() const -> const Eigen::VectorXd&
{
    throw GelexException("Result::allele_freq is legacy and will be removed");
}

}  // namespace gelex
