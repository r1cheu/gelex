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

#ifndef GELEX_ENGINE_PREDICT_H_
#define GELEX_ENGINE_PREDICT_H_

#include <filesystem>
#include <optional>

#include <Eigen/Core>

#include "gelex/data/dataframe/dataframe.h"
#include "gelex/infra/logging/predict_event.h"
#include "gelex/io/predict/input_reader.h"
#include "gelex/predict/snp_alignment.h"
#include "gelex/predict/types.h"

namespace gelex
{

class PredictEngine
{
   public:
    struct Config
    {
        std::string bfile_prefix;
        std::string gfile_prefix;

        std::optional<std::filesystem::path> qcovar_path;
        std::optional<std::filesystem::path> dcovar_path;
        std::filesystem::path output_path;
    };

    struct PredictParams
    {
        dataframe::DataFrame<std::string> snp_effects;
        Eigen::VectorXd add_effects;
        std::optional<Eigen::VectorXd> dom_effects;
        predict::Coefficients coefficients;
        predict::SbinData sbin;
    };

    struct PredictData
    {
        dataframe::DataFrame<std::string> fam_df;
        dataframe::DataFrame<std::string> bim_df;
        Eigen::MatrixXd covariates;
    };

    explicit PredictEngine(Config config);
    auto run(const PredictObserver& observer = {}) -> void;

   private:
    static auto load_sbin(const std::filesystem::path& path)
        -> predict::SbinData;

    auto load_params() const -> PredictParams;
    auto load_data(const PredictParams& params) const -> PredictData;
    auto align(
        const PredictParams& params,
        const PredictData& data,
        const PredictObserver& observer) const -> predict::SnpAlignment;
    auto select(
        const PredictData& data,
        const predict::SnpAlignment& alignment,
        bool has_dom) const -> predict::GenotypeData;
    Config config_;
};

}  // namespace gelex

#endif  // GELEX_ENGINE_PREDICT_H_
