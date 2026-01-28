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

#include <format>
#include <fstream>
#include <ostream>

#include "../src/data/parser.h"
#include "gelex/exception.h"

namespace gelex
{
namespace
{

void write_prediction_impl(
    std::ostream& stream,
    double total_prediction,
    const Eigen::Ref<const Eigen::RowVectorXd>& covar_pred)
{
    stream << std::format("\t{:.6f}", total_prediction);
    const Eigen::Index covar_end = covar_pred.cols();

    for (Eigen::Index j = 0; j < covar_end; ++j)
    {
        stream << std::format("\t{:.6f}", covar_pred(j));
    }
}

}  // namespace

PredictWriter::PredictWriter(
    const std::filesystem::path& output_path,
    bool iid_only)
    : output_path_(output_path), iid_only_(iid_only)
{
    if (output_path_.empty())
    {
        throw InvalidInputException("Output path must be provided");
    }
}

void PredictWriter::write_header(
    std::ostream& stream,
    std::span<const std::string> covar_names,
    bool has_dom) const
{
    if (!iid_only_)
    {
        stream << "FID\t";
    }
    stream << "IID\tprediction";

    for (const auto& name : covar_names)
    {
        stream << "\t" << name;
    }

    if (!has_dom)
    {
        stream << "\tadditive\n";
        return;
    }
    stream << "\tadditive\tdominant\n";
}

void PredictWriter::write_prediction(
    std::ostream& stream,
    double total_prediction,
    const Eigen::Ref<const Eigen::RowVectorXd>& covar_pred,
    double add_pred,
    double dom_pred)
{
    write_prediction_impl(stream, total_prediction, covar_pred);
    stream << std::format("\t{:.6f}\t{:.6f}\n", add_pred, dom_pred);
}

void PredictWriter::write_prediction(
    std::ostream& stream,
    double total_prediction,
    const Eigen::Ref<const Eigen::RowVectorXd>& covar_pred,
    double add_pred)
{
    write_prediction_impl(stream, total_prediction, covar_pred);
    stream << std::format("\t{:.6f}\n", add_pred);
}

void PredictWriter::write_id(std::ostream& stream, std::string_view sample_id)
    const
{
    if (iid_only_)
    {
        stream << sample_id;
    }
    else
    {
        auto pos = sample_id.find('_');
        if (pos != std::string_view::npos)
        {
            stream << sample_id.substr(0, pos) << '\t'
                   << sample_id.substr(pos + 1);
        }
    }
}

void PredictWriter::write(
    const Eigen::Ref<Eigen::VectorXd>& predictions,
    std::span<const std::string> sample_ids,
    const Eigen::Ref<Eigen::VectorXd>& add_pred,
    const Eigen::Ref<Eigen::VectorXd>& dom_pred,
    const Eigen::Ref<Eigen::MatrixXd>& covar_pred,
    std::span<const std::string> covar_names)
{
    if (sample_ids.size() != static_cast<size_t>(predictions.size()))
    {
        throw InvalidInputException(
            std::format(
                "Dimension mismatch: {} FIDs but {} predictions",
                sample_ids.size(),
                predictions.size()));
    }
    bool has_dom = dom_pred.size() > 0;

    auto stream = detail::open_file<std::ofstream>(output_path_, std::ios::out);
    write_header(stream, covar_names, has_dom);
    const Eigen::Index n_samples = predictions.size();

    if (has_dom)
    {
        for (Eigen::Index i = 0; i < n_samples; ++i)
        {
            write_id(stream, sample_ids[i]);
            write_prediction(
                stream,
                predictions[i],
                covar_pred.row(i),
                add_pred[i],
                dom_pred[i]);
        }
    }
    else
    {
        for (Eigen::Index i = 0; i < n_samples; ++i)
        {
            write_id(stream, sample_ids[i]);
            write_prediction(
                stream, predictions[i], covar_pred.row(i), add_pred[i]);
        }
    }
}

}  // namespace gelex
