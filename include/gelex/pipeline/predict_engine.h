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

#ifndef GELEX_PREDICT_PREDICT_ENGINE_H
#define GELEX_PREDICT_PREDICT_ENGINE_H

#include <filesystem>
#include <optional>

#include <Eigen/Core>

#include "gelex/data/dataframe/dataframe.h"
#include "gelex/infra/logging/predict_event.h"
#include "gelex/predict/predict_types.h"
#include "gelex/predict/reader.h"
#include "gelex/predict/snp_alignment.h"

namespace gelex
{

class PredictEngine
{
   public:
    struct Config
    {
        std::string bfile_prefix;
        std::string gfile_prefix;
        std::filesystem::path bed_path;

        std::optional<std::filesystem::path> qcovar_path;
        std::optional<std::filesystem::path> dcovar_path;
        std::filesystem::path output_path;
    };

    struct PredictParams
    {
        df::DataFrame<std::string> snp_effects;
        Eigen::VectorXd add_effects;
        std::optional<Eigen::VectorXd> dom_effects;
        Coefficients coefficients;
        SbinData sbin;
    };

    struct PredictData
    {
        df::DataFrame<std::string> fam_df;
        df::DataFrame<std::string> bim_df;
        Eigen::MatrixXd covariates;
    };

    explicit PredictEngine(Config config);
    auto run(const PredictObserver& observer = {}) -> void;

   private:
    static auto load_sbin(const std::filesystem::path& path) -> SbinData;

    auto load_params() const -> PredictParams;
    auto load_data(const PredictParams& params) const -> PredictData;
    auto align(
        const PredictParams& params,
        const PredictData& data,
        const PredictObserver& observer) const -> SnpAlignment;
    auto select(
        const PredictData& data,
        const SnpAlignment& alignment,
        bool has_dom) const -> GenotypeData;
    Config config_;
};

}  // namespace gelex

#endif  // GELEX_PREDICT_PREDICT_ENGINE_H
