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
#include <memory>
#include <span>
#include <string>
#include <string_view>

#include <Eigen/Core>

namespace gelex::detail
{
class TextWriter;
}

namespace gelex
{

class PredictWriter
{
   public:
    explicit PredictWriter(const std::filesystem::path& output_path);
    ~PredictWriter();

    auto write(
        const Eigen::Ref<const Eigen::VectorXd>& predictions,
        std::span<const std::string> sample_ids,
        const Eigen::Ref<const Eigen::VectorXd>& add_pred,
        const Eigen::Ref<const Eigen::VectorXd>& dom_pred,
        const Eigen::Ref<const Eigen::MatrixXd>& covar_pred,
        std::span<const std::string> covar_names) -> void;

   private:
    auto write_header(std::span<const std::string> covar_names, bool has_dom)
        -> void;

    auto write_prediction(
        double total_prediction,
        const Eigen::Ref<const Eigen::RowVectorXd>& covar_pred,
        double add_pred,
        bool has_dom,
        double dom_pred) -> void;

    auto write_id(std::string_view sample_id) -> void;

    std::unique_ptr<detail::TextWriter> writer_;
    std::string row_buf_;
};

}  // namespace gelex

#endif  // GELEX_PREDICT_PREDICT_WRITER_H
