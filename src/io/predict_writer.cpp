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

#include "gelex/io/predict_writer.h"

#include <fmt/format.h>
#include <Eigen/Core>
#include <cstddef>
#include <filesystem>
#include <memory>
#include <ranges>
#include <string>

#include "gelex/data/sample_id.h"
#include "gelex/exception.h"
#include "gelex/io/detail/text_writer.h"
#include "gelex/predict/types.h"
#include "gelex/types/genetic_mode.h"

namespace gelex
{

PredictWriter::PredictWriter(const std::filesystem::path& output_path)
{
    if (output_path.empty())
    {
        throw GelexException("Output path must be provided");
    }
    writer_ = std::make_unique<detail::TextWriter>(output_path);
}

PredictWriter::~PredictWriter() = default;

auto PredictWriter::write(const PredictResult& result) -> void
{
    const auto n_samples = static_cast<Eigen::Index>(result.sample_ids.size());
    if (n_samples != result.predictions.size())
    {
        throw GelexException(
            fmt::format(
                "Dimension mismatch: {} sample IDs but {} predictions",
                result.sample_ids.size(),
                result.predictions.size()));
    }

    std::string header = "FID\tIID\tprediction";
    for (const auto& name : result.covar_names)
    {
        header += '\t';
        header += name;
    }
    for (const auto& mode : std::views::keys(result.snp_components))
    {
        header += fmt::format("\t{}", mode);
    }
    writer_->write(header);

    for (Eigen::Index i = 0; i < n_samples; ++i)
    {
        row_buf_.clear();
        auto [fid, iid]
            = split_sample_id(result.sample_ids[static_cast<size_t>(i)]);
        row_buf_ += fmt::format("{}\t{}", fid, iid);
        row_buf_ += fmt::format("\t{:.6f}", result.predictions[i]);

        for (Eigen::Index j = 0; j < result.covar_predictions.cols(); ++j)
        {
            row_buf_ += fmt::format("\t{:.6f}", result.covar_predictions(i, j));
        }
        for (const auto& component : std::views::values(result.snp_components))
        {
            row_buf_ += fmt::format("\t{:.6f}", component[i]);
        }
        writer_->write(row_buf_);
    }
}

}  // namespace gelex
