#include "predict/predict_writer.h"

#include <format>
#include <fstream>
#include <ostream>
#include <ranges>

#include "../src/data/parser.h"
#include "gelex/exception.h"

namespace gelex
{

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
    std::span<const std::string> covar_names) const
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
    stream << "\tadditive\tdominant\n";
}
void PredictWriter::write_prediction_with_dom(
    std::ostream& stream,
    double total_prediction,
    const Eigen::Ref<const Eigen::RowVectorXd>& covar_pred,
    double add_pred,
    double dom_pred)
{
    stream << std::format("\t{:.6f}", total_prediction);
    const Eigen::Index covar_end = covar_pred.cols();

    for (Eigen::Index j = 0; j < covar_end; ++j)
    {
        stream << std::format("\t{:.6f}", covar_pred(j));
    }

    stream << std::format("\t{:.6f}\t{:.6f}\n", add_pred, dom_pred);
}

void PredictWriter::write_prediction_no_dom(
    std::ostream& stream,
    double total_prediction,
    const Eigen::Ref<const Eigen::RowVectorXd>& covar_pred,
    double add_pred)
{
    stream << std::format("\t{:.6f}", total_prediction);
    const Eigen::Index covar_end = covar_pred.cols();

    for (Eigen::Index j = 0; j < covar_end; ++j)
    {
        stream << std::format("\t{:.6f}", covar_pred(j));
    }

    stream << std::format("\t{:.6f}\n", add_pred);
}

void PredictWriter::write_id(std::ostream& stream, std::string_view sample_id)
{
    auto parts = sample_id | std::views::split('_') | std::views::take(2);
    auto it = parts.begin();

    if (iid_only_)
    {
        if (it != parts.end())
        {
            stream << std::string_view(*it);
        }
    }
    else
    {
        if (it != parts.end())
        {
            stream << std::string_view(*it);
            ++it;
        }
        if (it != parts.end())
        {
            stream << '\t' << std::string_view(*it);
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

    auto stream = detail::open_file<std::ofstream>(output_path_, std::ios::out);
    write_header(stream, covar_names);
    const Eigen::Index n_samples = predictions.size();

    if (dom_pred.size() > 0)
    {
        for (Eigen::Index i = 0; i < n_samples; ++i)
        {
            write_id(stream, sample_ids[i]);
            write_prediction_with_dom(
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
            write_prediction_no_dom(
                stream, predictions[i], covar_pred.row(i), add_pred[i]);
        }
    }
}

}  // namespace gelex
