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

#include "gelex/algo/infer/mcmc/records.h"

#include <string>
#include <string_view>
#include <type_traits>
#include <utility>
#include <variant>

#include <Eigen/Core>

#include "gelex/bayes/state.h"
#include "gelex/exception.h"
#include "gelex/infra/field_flag.h"

namespace gelex::mcmc
{

Records::CategoricalRecord::CategoricalRecord(
    Eigen::Index n_items,
    Eigen::Index n_categories)
    : draws(n_items, n_categories)
{
}

Records::Records(Records&& other) noexcept
    : records_(std::move(other.records_)),
      indices_(std::move(other.indices_)),
      category_counts_(std::move(other.category_counts_))
{
}

auto Records::operator=(Records&& other) noexcept -> Records&
{
    records_ = std::move(other.records_);
    indices_ = std::move(other.indices_);
    category_counts_ = std::move(other.category_counts_);
    return *this;
}

auto Records::store(BayesState& state) -> void
{
    state.visit(*this);
}

auto Records::result(std::string_view path) -> RecordResult
{
    auto key = std::string{path};
    const auto it = indices_.find(key);
    if (it == indices_.end())
    {
        throw GelexException("Records: missing record " + key);
    }
    auto& record = records_[it->second];
    return std::visit(
        []<typename T>(T& value) -> RecordResult
        {
            if constexpr (
                std::is_same_v<T, ScalarRecord>
                || std::is_same_v<T, VectorRecord>)
            {
                return value.draws.stats();
            }
            else
            {
                return std::move(value.draws).take_probabilities();
            }
        },
        record);
}

template <typename RecordType, typename Value>
auto Records::store_record(std::string_view name, Value&& value) -> void
{
    auto field_key = field_path(name);

    auto it = indices_.find(field_key);
    if (it == indices_.end())
    {
        const auto index = records_.size();
        records_.emplace_back(RecordType{});
        indices_.emplace(std::move(field_key), index);
        std::get<RecordType>(records_[index])
            .draws.store(std::forward<Value>(value));
        return;
    }

    auto& record = records_[it->second];
    if (!std::holds_alternative<RecordType>(record))
    {
        throw GelexException("Records: record kind changed for " + field_key);
    }
    std::get<RecordType>(record).draws.store(std::forward<Value>(value));
}

auto Records::on(
    std::string_view name,
    Eigen::Ref<Eigen::VectorXf> value,
    FieldFlag flags) -> void
{
    if (!has(flags, FieldFlag::trace))
    {
        return;
    }
    static_cast<void>(value);
    throw GelexException(
        "Records: VectorXf trace is not implemented for " + field_path(name));
}

auto Records::on(
    std::string_view name,
    Eigen::Ref<Eigen::VectorXd> value,
    FieldFlag flags) -> void
{
    if (!has(flags, FieldFlag::trace))
    {
        return;
    }

    if (name == "proportion")
    {
        category_counts_[std::string{path()}] = value.size();
    }
    store_record<VectorRecord>(name, value);
}

auto Records::on(
    std::string_view name,
    Eigen::Ref<Eigen::VectorXi> value,
    FieldFlag flags) -> void
{
    if (!has(flags, FieldFlag::trace))
    {
        return;
    }
    if (name != "assignment")
    {
        throw GelexException(
            "Records: traced integer vector " + std::string{path()});
    }

    auto scope_key = std::string{path()};
    auto category_it = category_counts_.find(scope_key);
    if (category_it == category_counts_.end())
    {
        throw GelexException(
            "Records: assignment without proportion " + scope_key);
    }

    auto field_key = field_path(name);

    auto it = indices_.find(field_key);
    if (it == indices_.end())
    {
        const auto index = records_.size();
        records_.emplace_back(
            CategoricalRecord{value.size(), category_it->second});
        indices_.emplace(std::move(field_key), index);
        std::get<CategoricalRecord>(records_[index]).draws.store(value);
        return;
    }

    auto& record = records_[it->second];
    if (!std::holds_alternative<CategoricalRecord>(record))
    {
        throw GelexException("Records: record kind changed for " + field_key);
    }
    std::get<CategoricalRecord>(record).draws.store(value);
}

auto Records::on(std::string_view name, double& value, FieldFlag flags) -> void
{
    if (!has(flags, FieldFlag::trace))
    {
        return;
    }

    store_record<ScalarRecord>(name, value);
}

auto Records::on(std::string_view name, int& value, FieldFlag flags) -> void
{
    if (!has(flags, FieldFlag::trace))
    {
        return;
    }

    static_cast<void>(value);
    throw GelexException(
        "Records: int trace is not implemented for " + field_path(name));
}

}  // namespace gelex::mcmc
