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
        std::span<const std::string> covar_names,
        bool has_dom) const;

    static void write_prediction(
        std::ostream& stream,
        double total_prediction,
        const Eigen::Ref<const Eigen::RowVectorXd>& covar_pred,
        double add_pred,
        double dom_pred);

    static void write_prediction(
        std::ostream& stream,
        double total_prediction,
        const Eigen::Ref<const Eigen::RowVectorXd>& covar_pred,
        double add_pred);

    void write_id(std::ostream& stream, std::string_view sample_id) const;

    std::filesystem::path output_path_;
    bool iid_only_ = false;
};

}  // namespace gelex

#endif  // GELEX_PREDICT_PREDICT_WRITER_H
