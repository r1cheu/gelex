#ifndef GELEX_PREDICT_PREDICT_WRITER_H
#define GELEX_PREDICT_PREDICT_WRITER_H

#include <filesystem>
#include <string>
#include <vector>

#include <Eigen/Core>

namespace gelex
{

class PredictWriter
{
   public:
    explicit PredictWriter(const std::filesystem::path& output_path);

    void write(
        const Eigen::VectorXd& predictions,
        const std::vector<std::string>& fids,
        const std::vector<std::string>& iids,
        const Eigen::VectorXd& add_pred,
        const Eigen::VectorXd& dom_pred,
        const Eigen::MatrixXd& covar_pred,
        const std::vector<std::string>& covar_names);

   private:
    static void write_header(
        std::ostream& stream,
        const std::vector<std::string>& covar_names);

    static void write_sample(
        std::ostream& stream,
        Eigen::Index sample_idx,
        const std::string& fid,
        const std::string& iid,
        double total_prediction,
        const Eigen::MatrixXd& covar_pred,
        double add_pred,
        double dom_pred);

    std::filesystem::path output_path_;
};

}  // namespace gelex

#endif  // GELEX_PREDICT_PREDICT_WRITER_H
