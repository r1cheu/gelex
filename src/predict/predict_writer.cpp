#include "predict/predict_writer.h"

#include <format>
#include <fstream>
#include <ostream>

#include "../src/data/parser.h"
#include "gelex/exception.h"

namespace gelex
{

PredictWriter::PredictWriter(const std::filesystem::path& output_path)
    : output_path_(output_path)
{
    if (output_path_.empty())
    {
        throw InvalidInputException("Output path must be provided");
    }
}

void PredictWriter::write_header(
    std::ostream& stream,
    const std::vector<std::string>& covar_names)
{
    stream << "FID\tIID\tprediction";
    for (const auto& name : covar_names)
    {
        stream << "\t" << name;
    }
    stream << "\tadditive\tdominant\n";
}

void PredictWriter::write_sample(
    std::ostream& stream,
    Eigen::Index sample_idx,
    const std::string& fid,
    const std::string& iid,
    double total_prediction,
    const Eigen::MatrixXd& covar_pred,
    double add_pred,
    double dom_pred)
{
    stream << fid << "\t" << iid;
    stream << std::format("\t{:.6f}", total_prediction);

    const Eigen::Index covar_end = covar_pred.cols();

    for (Eigen::Index j = 0; j < covar_end; ++j)
    {
        stream << std::format("\t{:.6f}", covar_pred(sample_idx, j));
    }

    stream << std::format("\t{:.6f}\t{:.6f}", add_pred, dom_pred);
    stream << "\n";
}

void PredictWriter::write(
    const Eigen::VectorXd& predictions,
    const std::vector<std::string>& fids,
    const std::vector<std::string>& iids,
    const Eigen::VectorXd& add_pred,
    const Eigen::VectorXd& dom_pred,
    const Eigen::MatrixXd& covar_pred,
    const std::vector<std::string>& covar_names)
{
    if (fids.size() != static_cast<size_t>(predictions.size()))
    {
        throw InvalidInputException(
            std::format(
                "Dimension mismatch: {} FIDs but {} predictions",
                fids.size(),
                predictions.size()));
    }
    if (iids.size() != static_cast<size_t>(predictions.size()))
    {
        throw InvalidInputException(
            std::format(
                "Dimension mismatch: {} IIDs but {} predictions",
                iids.size(),
                predictions.size()));
    }

    auto stream = detail::open_file<std::ofstream>(output_path_, std::ios::out);

    write_header(stream, covar_names);

    const Eigen::Index n_samples = predictions.size();

    for (Eigen::Index i = 0; i < n_samples; ++i)
    {
        write_sample(
            stream,
            i,
            fids[i],
            iids[i],
            predictions[i],
            covar_pred,
            add_pred[i],
            dom_pred[i]);
    }
}

}  // namespace gelex
