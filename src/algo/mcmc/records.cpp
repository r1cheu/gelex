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

#include "gelex/algo/mcmc/records.h"

#include <Eigen/Core>
#include <memory>
#include <ranges>
#include <string>
#include <string_view>
#include <type_traits>
#include <utility>
#include <variant>

#include "gelex/bayes/legacy_state.h"
#include "gelex/bayes/model.h"
#include "gelex/exception.h"
#include "gelex/infra/field_flag.h"
#include "gelex/io/binary_writer.h"

namespace gelex
{

Records::Records(Eigen::Index n_draws, std::string_view draws_path)
    : n_draws_(n_draws)
{
    if (!draws_path.empty())
    {
        if (n_draws_ <= 0)
        {
            throw GelexException("Records: n_draws must be positive");
        }
        writer_ = std::make_unique<BinaryWriter>(draws_path);
    }
}

Records::Records(Records&& other) noexcept
    : records_(std::move(other.records_)),
      paths_(std::move(other.paths_)),
      names_(std::move(other.names_)),
      indices_(std::move(other.indices_)),
      writer_(std::move(other.writer_)),
      n_draws_(other.n_draws_)
{
}

auto Records::operator=(Records&& other) noexcept -> Records&
{
    records_ = std::move(other.records_);
    paths_ = std::move(other.paths_);
    names_ = std::move(other.names_);
    indices_ = std::move(other.indices_);
    writer_ = std::move(other.writer_);
    n_draws_ = other.n_draws_;
    return *this;
}

Records::~Records() = default;

auto Records::store(const BayesModel& model, LegacyBayesState& state) -> void
{
    state.visit(*this);
    model.visit(*this);
}

auto Records::take_results() && -> std::vector<RecordEntry>
{
    std::vector<RecordEntry> output;
    output.reserve(paths_.size());
    for (auto [i, path] : std::views::enumerate(paths_))
    {
        auto value = result(path);
        output.push_back(
            RecordEntry{
                std::move(path),
                std::move(value),
                std::move(names_[static_cast<std::size_t>(i)])});
    }

    records_.clear();
    if (writer_)
    {
        writer_->close();
        writer_.reset();
    }

    paths_.clear();
    names_.clear();
    indices_.clear();
    return output;
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
                std::is_same_v<T, detail::ScalarRecord>
                || std::is_same_v<T, detail::VectorRecord>)
            {
                return value.result();
            }
            else
            {
                return std::move(value).result();
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
        if constexpr (std::is_same_v<RecordType, detail::VectorRecord>)
        {
            records_.emplace_back(RecordType{*this, field_key, value});
        }
        else
        {
            records_.emplace_back(RecordType{*this, field_key});
        }
        paths_.push_back(field_key);
        names_.emplace_back(std::nullopt);
        auto& stored = std::get<RecordType>(records_[index]);
        indices_.emplace(std::move(field_key), index);
        stored.store(std::forward<Value>(value));
        return;
    }

    auto& record = records_[it->second];
    if (!std::holds_alternative<RecordType>(record))
    {
        throw GelexException("Records: record kind changed for " + field_key);
    }
    auto& stored = std::get<RecordType>(record);
    stored.store(std::forward<Value>(value));
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

    store_record<detail::VectorRecord>(name, value);
}

auto Records::on(
    std::string_view name,
    CategoricalVector& value,
    FieldFlag flags) -> void
{
    if (!has(flags, FieldFlag::trace))
    {
        return;
    }
    if (value.category_count() <= 0)
    {
        throw GelexException(
            "Records: category count must be positive for " + field_path(name));
    }
    auto field_key = field_path(name);

    auto it = indices_.find(field_key);
    if (it == indices_.end())
    {
        const auto index = records_.size();
        records_.emplace_back(
            detail::CategoricalRecord{*this, field_key, value});
        paths_.push_back(field_key);
        names_.emplace_back(std::nullopt);
        auto& stored = std::get<detail::CategoricalRecord>(records_[index]);
        indices_.emplace(std::move(field_key), index);
        stored.store(value);
        return;
    }

    auto& record = records_[it->second];
    if (!std::holds_alternative<detail::CategoricalRecord>(record))
    {
        throw GelexException("Records: record kind changed for " + field_key);
    }
    auto& stored = std::get<detail::CategoricalRecord>(record);
    stored.store(value);
}

auto Records::on(std::string_view name, double& value, FieldFlag flags) -> void
{
    if (!has(flags, FieldFlag::trace))
    {
        return;
    }

    store_record<detail::ScalarRecord>(name, value);
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

auto Records::on(
    std::string_view name,
    std::span<const std::string> value,
    FieldFlag flags) -> void
{
    if (!has(flags, FieldFlag::report))
    {
        return;
    }
    if (value.empty())
    {
        return;
    }

    constexpr std::string_view suffix{"_names"};
    if (!name.ends_with(suffix))
    {
        return;
    }

    const auto field_key
        = field_path(name.substr(0, name.size() - suffix.size()));
    const auto it = indices_.find(field_key);
    if (it == indices_.end())
    {
        return;
    }
    if (names_[it->second])
    {
        return;
    }
    names_[it->second].emplace(value.begin(), value.end());
}

auto Records::on(std::string_view name, std::string_view value, FieldFlag flags)
    -> void
{
    if (!has(flags, FieldFlag::report))
    {
        return;
    }
    if (value.empty())
    {
        return;
    }

    constexpr std::string_view suffix{"_name"};
    if (!name.ends_with(suffix))
    {
        return;
    }

    const auto field_key
        = field_path(name.substr(0, name.size() - suffix.size()));
    const auto it = indices_.find(field_key);
    if (it == indices_.end())
    {
        return;
    }
    if (names_[it->second])
    {
        return;
    }
    names_[it->second].emplace(1, std::string{value});
}

}  // namespace gelex
