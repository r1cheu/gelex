#ifndef GELEX_PREDICT_PREDICT_WRITER_H
#define GELEX_PREDICT_PREDICT_WRITER_H

#include <filesystem>
#include <span>
#include <string>

#include <Eigen/Core>

namespace gelex
{

class PredictWriter
{
   public:
    PredictWriter(const std::filesystem::path& output_path, bool iid_only);

    void write(
        const Eigen::Ref<Eigen::VectorXd>& predictions,
        std::span<const std::string> sample_ids,
        const Eigen::Ref<Eigen::VectorXd>& add_pred,
        const Eigen::Ref<Eigen::VectorXd>& dom_pred,
        const Eigen::Ref<Eigen::MatrixXd>& covar_pred,
        std::span<const std::string> covar_names);

   private:
    void write_header(
        std::ostream& stream,
        std::span<const std::string> covar_names) const;

    void write_prediction_with_dom(
        std::ostream& stream,
        double total_prediction,
        const Eigen::Ref<const Eigen::RowVectorXd>& covar_pred,
        double add_pred,
        double dom_pred);

    void write_prediction_no_dom(
        std::ostream& stream,
        double total_prediction,
        const Eigen::Ref<const Eigen::RowVectorXd>& covar_pred,
        double add_pred);

    // split sample by "_" use ranges
    void write_id(std::ostream& stream, std::string_view sample_id);

    std::filesystem::path output_path_;
    bool iid_only_ = false;
};

}  // namespace gelex

#endif  // GELEX_PREDICT_PREDICT_WRITER_H
