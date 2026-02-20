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

#include "pipeline/predict/predict_writer.h"

#include <cstddef>
#include <format>
#include <memory>

#include "gelex/data/frame/dataframe_policy.h"
#include "gelex/exception.h"
#include "gelex/io/text_writer.h"

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
    std::span<const std::string> covar_names,
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

auto PredictWriter::write_prediction(
    double total_prediction,
    const Eigen::Ref<const Eigen::RowVectorXd>& covar_pred,
    double add_pred,
    bool has_dom,
    double dom_pred) -> void
{
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
}

auto PredictWriter::write_id(std::string_view sample_id) -> void
{
    auto [fid, iid] = split_sample_id(sample_id);
    row_buf_ += std::format("{}\t{}", fid, iid);
}

auto PredictWriter::write(
    const Eigen::Ref<const Eigen::VectorXd>& predictions,
    std::span<const std::string> sample_ids,
    const Eigen::Ref<const Eigen::VectorXd>& add_pred,
    const Eigen::Ref<const Eigen::VectorXd>& dom_pred,
    const Eigen::Ref<const Eigen::MatrixXd>& covar_pred,
    std::span<const std::string> covar_names) -> void
{
    if (sample_ids.size() != static_cast<std::size_t>(predictions.size()))
    {
        throw InvalidInputException(
            std::format(
                "Dimension mismatch: {} FIDs but {} predictions",
                sample_ids.size(),
                predictions.size()));
    }

    const bool has_dom = dom_pred.size() > 0;
    write_header(covar_names, has_dom);

    const Eigen::Index n_samples = predictions.size();
    for (Eigen::Index i = 0; i < n_samples; ++i)
    {
        row_buf_.clear();
        write_id(sample_ids[i]);
        const double dom_value = has_dom ? dom_pred[i] : 0.0;
        write_prediction(
            predictions[i], covar_pred.row(i), add_pred[i], has_dom, dom_value);
        writer_->write(row_buf_);
    }
}

}  // namespace gelex
