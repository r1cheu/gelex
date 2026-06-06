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

#include <utility>

#include "gelex/algo/mcmc/records.h"
#include "gelex/io/binary_writer.h"

namespace gelex::mcmc::detail
{

ScalarRecord::ScalarRecord(Records& owner, std::string_view path)
{
    if (owner.writer_)
    {
        draw_handle_
            = owner.writer_->reserve<double>(path, 1, owner.n_draws_).index;
    }
}

auto ScalarRecord::store(Records& owner, double value) -> void
{
    draws_.store(value);
    if (owner.writer_)
    {
        owner.writer_->write(io::SectionHandle<double>{*draw_handle_}, value);
    }
}

auto ScalarRecord::result() const -> stats::RunningStatsResult
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
        draw_handle_
            = owner.writer_->reserve<double>(path, value.size(), owner.n_draws_)
                  .index;
    }
}

auto VectorRecord::store(
    Records& owner,
    const Eigen::Ref<const Eigen::VectorXd>& value) -> void
{
    draws_.store(value);
    if (owner.writer_)
    {
        owner.writer_->write(io::SectionHandle<double>{*draw_handle_}, value);
    }
}

auto VectorRecord::result() const -> stats::RunningStatsResult
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
        draw_handle_
            = owner.writer_->reserve<int>(path, value.size(), owner.n_draws_)
                  .index;
    }
}

auto CategoricalRecord::store(Records& owner, const CategoricalVector& value)
    -> void
{
    draws_.store(value);
    if (owner.writer_)
    {
        owner.writer_->write(io::SectionHandle<int>{*draw_handle_}, value);
    }
}

auto CategoricalRecord::result() && -> stats::CategoryProbResult
{
    return std::move(draws_).take_probabilities();
}

}  // namespace gelex::mcmc::detail
