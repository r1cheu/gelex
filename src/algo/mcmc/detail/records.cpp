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

#include "gelex/algo/mcmc/detail/records.h"

#include <cstddef>
#include <cstdint>
#include <span>
#include <utility>

#include "gelex/algo/mcmc/records.h"
#include "gelex/io/binary_writer.h"

namespace gelex::detail
{

ScalarRecord::ScalarRecord(Records& owner, std::string_view path)
{
    if (owner.writer_)
    {
        draw_payload_ = owner.writer_->reserve<double>(
            path, BinaryShape{1, static_cast<std::uint64_t>(owner.n_draws_)});
    }
}

auto ScalarRecord::store(double value) -> void
{
    draws_.store(value);
    if (draw_payload_)
    {
        draw_payload_->append(value);
    }
}

auto ScalarRecord::result() const -> RunningStatsResult
{
    return draws_.stats();
}

VectorRecord::VectorRecord(
    Records& owner,
    std::string_view path,
    const Eigen::Ref<const Eigen::VectorXd>& value)
{
    if (owner.writer_)
    {
        draw_payload_ = owner.writer_->reserve<double>(
            path,
            BinaryShape{
                static_cast<std::uint64_t>(value.size()),
                static_cast<std::uint64_t>(owner.n_draws_)});
    }
}

auto VectorRecord::store(const Eigen::Ref<const Eigen::VectorXd>& value) -> void
{
    draws_.store(value);
    if (draw_payload_)
    {
        draw_payload_->append(
            std::span<const double>{
                value.data(), static_cast<std::size_t>(value.size())});
    }
}

auto VectorRecord::result() const -> RunningStatsResult
{
    return draws_.stats();
}

CategoricalRecord::CategoricalRecord(
    Records& owner,
    std::string_view path,
    const CategoricalVector& value)
    : draws_(value.size(), value.category_count())
{
    if (owner.writer_)
    {
        draw_payload_ = owner.writer_->reserve<int>(
            path,
            BinaryShape{
                static_cast<std::uint64_t>(value.size()),
                static_cast<std::uint64_t>(owner.n_draws_)});
    }
}

auto CategoricalRecord::store(const CategoricalVector& value) -> void
{
    draws_.store(value);
    if (draw_payload_)
    {
        draw_payload_->append(
            std::span<const int>{
                value.data(), static_cast<std::size_t>(value.size())});
    }
}

auto CategoricalRecord::result() && -> CategoryProbResult
{
    return std::move(draws_).take_probabilities();
}

}  // namespace gelex::detail
