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

#include "../src/predict/genotype_aligner.h"

#include <Eigen/Core>

#include "../src/predict/snp_matcher.h"

namespace gelex
{

GenotypeAligner::GenotypeAligner(
    const std::filesystem::path& bed_path,
    const SnpEffects& snp_effects)
    : num_snp_effects_(static_cast<Eigen::Index>(snp_effects.size()))
{
    SnpMatcher matcher(snp_effects);
    match_plan_ = matcher.match(bed_path);
}

auto GenotypeAligner::align(Eigen::MatrixXd&& raw_genotype) const
    -> Eigen::MatrixXd
{
    Eigen::MatrixXd genotype = std::move(raw_genotype);
    const Eigen::Index num_samples = genotype.rows();

    Eigen::MatrixXd result(num_samples, num_snp_effects_);
    result.setZero();

#pragma omp parallel for schedule(guided) default(none) shared(result, genotype)
    for (Eigen::Index i = 0; i < static_cast<Eigen::Index>(match_plan_.size());
         ++i)
    {
        const MatchInfo& info = match_plan_[i];
        if (info.type == MatchType::skip || info.target_col < 0)
        {
            continue;
        }
        if (info.type == MatchType::reverse)
        {
            result.col(info.target_col) = 2.0 - genotype.col(i).array();
        }
        else
        {
            result.col(info.target_col) = genotype.col(i);
        }
    }
    return result;
}
}  // namespace gelex
