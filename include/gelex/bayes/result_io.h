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

#ifndef GELEX_BAYES_RESULT_IO_H_
#define GELEX_BAYES_RESULT_IO_H_

#include <fmt/format.h>
#include <string_view>

#include "gelex/bayes/detail/result_writer.h"
#include "gelex/bayes/result.h"
#include "gelex/genetic_mode.h"
#include "gelex/io/detail/text_writer.h"

namespace gelex
{

template <typename GeneticPrior>
auto write_params(
    const BayesResult<GeneticPrior>& result,
    std::string_view prefix) -> void
{
    detail::TextWriter writer{fmt::format("{}.params", prefix)};
    writer.write_header({"term", "index", "column_name", "mean", "stddev"});

    detail::write_parameter_rows(writer, result.fixed());
    for (const auto& random : result.random())
    {
        detail::write_parameter_rows(writer, random.coefficients());
    }
}

template <typename GeneticPrior>
auto write_summary(
    const BayesResult<GeneticPrior>& result,
    std::string_view prefix) -> void
{
    detail::TextWriter writer{fmt::format("{}.summary", prefix)};
    writer.write_header({"term", "index", "mean", "stddev"});

    for (const auto& random : result.random())
    {
        detail::write_summary_rows(writer, random.variance());
        detail::write_summary_rows(writer, random.explained_variance());
    }

    detail::write_genetic_summary_rows(writer, result.genetic());
    detail::write_summary_rows(writer, result.residual());

    const auto& variance_summary = result.variance_summary();
    const auto write_mode_summary = [&writer]<GeneticMode>(const auto& summary)
    { detail::write_summary_rows(writer, summary); };
    variance_summary.explained_variances().for_each(write_mode_summary);
    variance_summary.heritabilities().for_each(write_mode_summary);
    detail::write_summary_rows(
        writer, variance_summary.total_explained_variance());
    detail::write_summary_rows(writer, variance_summary.total_heritability());
}

}  // namespace gelex

#endif  // GELEX_BAYES_RESULT_IO_H_
