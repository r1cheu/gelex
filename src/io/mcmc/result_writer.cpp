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

#include "gelex/io/mcmc/result_writer.h"

#include <cstddef>
#include <filesystem>
#include <variant>

#include <fmt/format.h>
#include <Eigen/Core>

#include "gelex/algo/infer/mcmc/result.h"
#include "gelex/infra/stats/result.h"
#include "gelex/io/detail/text_writer.h"

namespace gelex
{

auto mcmc::write_result(
    const Result& result,
    const std::filesystem::path& prefix) -> void
{
    auto output_path = prefix;
    output_path += ".summary";

    io::detail::TextWriter writer(output_path);
    writer.write_header({"term", "mean", "stddev"});

    for (const auto& record : result.records())
    {
        if (!record.names)
        {
            continue;
        }

        if (!std::holds_alternative<stats::RunningStatsResult>(record.value))
        {
            continue;
        }

        const auto& stats = std::get<stats::RunningStatsResult>(record.value);
        const auto& names = *record.names;
        for (Eigen::Index i = 0; i < static_cast<Eigen::Index>(names.size());
             ++i)
        {
            writer.write(
                fmt::format(
                    "{}\t{}\t{}",
                    names[static_cast<std::size_t>(i)],
                    stats.mean(i),
                    stats.stddev(i)));
        }
    }
}

}  // namespace gelex
