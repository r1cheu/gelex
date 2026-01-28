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

#ifndef GELEX_PREDICT_GENOTYPE_ALIGNER_H
#define GELEX_PREDICT_GENOTYPE_ALIGNER_H

#include <Eigen/Core>

#include "../src/predict/snp_matcher.h"

namespace gelex
{
class GenotypeAligner
{
   public:
    GenotypeAligner(
        const std::filesystem::path& bed_path,
        const SnpEffects& snp_effects);

    auto align(Eigen::MatrixXd&& raw_genotype) const -> Eigen::MatrixXd;

   private:
    MatchPlan match_plan_;
    Eigen::Index num_snp_effects_ = 0;
};

}  // namespace gelex

#endif  // GELEX_PREDICT_GENOTYPE_ALIGNER_H
