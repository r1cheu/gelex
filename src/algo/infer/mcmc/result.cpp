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

#include "gelex/algo/infer/mcmc/records.h"
#include "gelex/bayes/model.h"
#include "gelex/exception.h"

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

}  // namespace gelex
