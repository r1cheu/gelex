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

#ifndef GELEX_PREDICT_SNP_PREDICTOR_H
#define GELEX_PREDICT_SNP_PREDICTOR_H

#include <Eigen/Core>

#include "gelex/types/snp_info.h"

namespace gelex
{

struct SnpComputeResult
{
    Eigen::VectorXd add;
    Eigen::VectorXd dom;
    [[nodiscard]] Eigen::VectorXd total() const
    {
        if (dom.size() > 0)
        {
            return add + dom;
        }
        return add;
    }
};

class SnpPredictor
{
   public:
    explicit SnpPredictor(const SnpEffects& effects);

    SnpComputeResult compute(
        const Eigen::Ref<const Eigen::MatrixXd>& genotype) const;

   private:
    SnpEffects effects_;
    void validate_dimensions(
        const Eigen::Ref<const Eigen::MatrixXd>& genotype) const;
};

}  // namespace gelex

#endif  // GELEX_PREDICT_SNP_PREDICTOR_H
