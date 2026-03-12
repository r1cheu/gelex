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

#ifndef GELEX_PREDICT_GENOTYPE_SELECTOR_H
#define GELEX_PREDICT_GENOTYPE_SELECTOR_H

#include <filesystem>
#include <vector>

#include <Eigen/Core>

#include "gelex/infra/logging/predict_event.h"
#include "gelex/types/snp_effects.h"

namespace gelex
{

class GenotypeSelector
{
   public:
    GenotypeSelector(
        const std::filesystem::path& bed_path,
        const SnpEffects& snp_effects,
        const PredictObserver& observer = {});

    auto select(Eigen::MatrixXd&& raw_genotype) const -> Eigen::MatrixXd;

   private:
    static constexpr auto normalize_allele(char allele) noexcept -> char
    {
        return static_cast<char>(allele & 0xDF);
    }

    auto build_column_map(const SnpIndex& bim) -> void;
    auto check_missing(const PredictObserver& observer) -> void;

    // -1 sentinel means missing SNP (mean-imputed in select())
    std::vector<Eigen::Index> column_map_;
    std::vector<Eigen::Index> missing_indices_;
    const SnpEffects* snp_effects_ = nullptr;
    Eigen::Index num_effect_snps_ = 0;
};

}  // namespace gelex

#endif  // GELEX_PREDICT_GENOTYPE_SELECTOR_H
