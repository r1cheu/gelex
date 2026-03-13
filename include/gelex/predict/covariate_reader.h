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

#ifndef GELEX_PREDICT_COVARIATE_READER_H
#define GELEX_PREDICT_COVARIATE_READER_H

#include <filesystem>
#include <string>
#include <vector>

#include <Eigen/Core>

namespace gelex
{

struct CovariateParams
{
    std::vector<std::string> term_names;
    Eigen::VectorXd coefficients;
};

class CovariateReader
{
   public:
    explicit CovariateReader(const std::filesystem::path& param_file_path)
        : params_(parse(param_file_path))
    {
    }

    auto params() const -> const CovariateParams& { return params_; }
    auto take_params() && -> CovariateParams { return std::move(params_); }

   private:
    static auto parse(const std::filesystem::path& path) -> CovariateParams;
    CovariateParams params_;
};

}  // namespace gelex

#endif  // GELEX_PREDICT_COVARIATE_READER_H
