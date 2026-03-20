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

#include "predict/predict_writer.h"

#include <cstddef>
#include <format>

#include "gelex/exception.h"
#include "gelex/io/text_writer.h"
#include "gelex/types/sample_id.h"

namespace gelex
{

PredictWriter::PredictWriter(const std::filesystem::path& output_path)
{
    if (output_path.empty())
    {
        throw InvalidInputException("Output path must be provided");
    }
    writer_ = std::make_unique<detail::TextWriter>(output_path);
}

PredictWriter::~PredictWriter() = default;

auto PredictWriter::write_header(
    const std::vector<std::string>& covar_names,
    bool has_dom) -> void
{
    std::string h = "FID\tIID\tprediction";

    for (const auto& name : covar_names)
    {
        h += '\t';
        h += name;
    }

    h += "\tadditive";
    if (has_dom)
    {
        h += "\tdominant";
    }

    writer_->write(h);
}

auto PredictWriter::write_row(
    std::string_view sample_id,
    double total_prediction,
    const Eigen::Ref<const Eigen::RowVectorXd>& covar_pred,
    double add_pred,
    bool has_dom,
    double dom_pred) -> void
{
    row_buf_.clear();
    auto [fid, iid] = split_sample_id(sample_id);
    row_buf_ += std::format("{}\t{}", fid, iid);
    row_buf_ += std::format("\t{:.6f}", total_prediction);

    for (Eigen::Index j = 0; j < covar_pred.cols(); ++j)
    {
        row_buf_ += std::format("\t{:.6f}", covar_pred(j));
    }

    row_buf_ += std::format("\t{:.6f}", add_pred);
    if (has_dom)
    {
        row_buf_ += std::format("\t{:.6f}", dom_pred);
    }

    writer_->write(row_buf_);
}

auto PredictWriter::write(const PredictResult& result) -> void
{
    const auto n_samples = static_cast<Eigen::Index>(result.sample_ids.size());
    if (n_samples != result.predictions.size())
    {
        throw InvalidInputException(
            std::format(
                "Dimension mismatch: {} sample IDs but {} predictions",
                result.sample_ids.size(),
                result.predictions.size()));
    }

    const bool has_dom = result.dom_predictions.has_value();
    write_header(result.covar_names, has_dom);

    for (Eigen::Index i = 0; i < n_samples; ++i)
    {
        const double dom_value = has_dom ? (*result.dom_predictions)[i] : 0.0;
        write_row(
            result.sample_ids[static_cast<size_t>(i)],
            result.predictions[i],
            result.covar_predictions.row(i),
            result.add_predictions[i],
            has_dom,
            dom_value);
    }
}

}  // namespace gelex
