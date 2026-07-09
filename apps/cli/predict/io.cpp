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

#include "io.h"

#include <Eigen/Core>
#include <cstddef>
#include <filesystem>
#include <fmt/format.h>
#include <ranges>
#include <string>

#include "gelex/data/sample_id.h"
#include "gelex/exception.h"
#include "gelex/io/detail/text_writer.h"

namespace cli
{

auto write_predictions(
    const std::filesystem::path& output_path,
    std::span<const std::string> sample_ids,
    const Eigen::Ref<const Eigen::VectorXd>& prediction,
    const CovariateResult& covar,
    const gelex::ModeMap<Eigen::VectorXd>& gebvs) -> void
{
    if (output_path.empty())
    {
        throw gelex::GelexException("Output path must be provided");
    }

    const auto n_samples = static_cast<Eigen::Index>(sample_ids.size());
    if (n_samples != prediction.size())
    {
        throw gelex::GelexException(
            fmt::format(
                "Dimension mismatch: {} sample IDs but {} predictions",
                sample_ids.size(),
                prediction.size()));
    }

    gelex::detail::TextWriter writer(output_path);

    std::string header = "FID\tIID\tprediction";
    for (const auto& name : covar.covar_names)
    {
        header += '\t';
        header += name;
    }
    for (const auto& mode : std::views::keys(gebvs))
    {
        header += fmt::format("\t{}", mode);
    }
    writer.write(header);

    std::string row_buf;
    for (Eigen::Index i = 0; i < n_samples; ++i)
    {
        row_buf.clear();
        auto [fid, iid]
            = gelex::split_sample_id(sample_ids[static_cast<std::size_t>(i)]);
        row_buf += fmt::format("{}\t{}", fid, iid);
        row_buf += fmt::format("\t{:.6f}", prediction[i]);

        for (Eigen::Index j = 0; j < covar.per_covariate.cols(); ++j)
        {
            row_buf += fmt::format("\t{:.6f}", covar.per_covariate(i, j));
        }
        for (const auto& component : std::views::values(gebvs))
        {
            row_buf += fmt::format("\t{:.6f}", component[i]);
        }
        writer.write(row_buf);
    }
}

}  // namespace cli
