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

#include "gelex/predict/genotype_selector.h"

#include <format>
#include <string>
#include <unordered_map>
#include <vector>

#include <Eigen/Core>

#include "gelex/data/reader/bim_reader.h"
#include "gelex/exception.h"
#include "gelex/infra/logging/notify.h"
#include "gelex/types/snp_effects.h"

namespace gelex
{

static constexpr double kMaxMissingFraction = 0.2;

GenotypeSelector::GenotypeSelector(
    const std::filesystem::path& bed_path,
    const SnpEffects& snp_effects,
    const PredictObserver& observer)
    : snp_effects_(&snp_effects)
{
    auto bim_path = bed_path;
    bim_path.replace_extension(".bim");

    detail::BimReader bim_reader(bim_path);
    build_column_map(bim_reader.info());
    check_missing(observer);
}

auto GenotypeSelector::build_column_map(const SnpIndex& predict_bim) -> void
{
    std::unordered_map<std::string, Eigen::Index> bim_id_map;
    bim_id_map.reserve(predict_bim.size());
    for (size_t i = 0; i < predict_bim.size(); ++i)
    {
        bim_id_map.emplace(predict_bim[i].id, static_cast<Eigen::Index>(i));
    }

    num_effect_snps_ = static_cast<Eigen::Index>(snp_effects_->size());
    column_map_.resize(static_cast<size_t>(num_effect_snps_));

    for (size_t i = 0; i < snp_effects_->size(); ++i)
    {
        const auto& effect_meta = (*snp_effects_)[i];
        auto it = bim_id_map.find(effect_meta.id);

        if (it == bim_id_map.end())
        {
            column_map_[i] = -1;
            missing_indices_.push_back(static_cast<Eigen::Index>(i));
            continue;
        }

        const auto bim_col = it->second;
        const auto& bim_meta = predict_bim[static_cast<size_t>(bim_col)];

        const char eff_a1 = normalize_allele(effect_meta.A1);
        const char eff_a2 = normalize_allele(effect_meta.A2);
        const char bim_a1 = normalize_allele(bim_meta.A1);
        const char bim_a2 = normalize_allele(bim_meta.A2);

        if (eff_a1 == bim_a1 && eff_a2 == bim_a2)
        {
            column_map_[i] = bim_col;
        }
        else
        {
            throw InvalidInputException(
                std::format(
                    "Allele mismatch for SNP '{}': effect file has A1={}, "
                    "A2={} "
                    "but prediction data has A1={}, A2={}.\n"
                    "Please pre-align with: plink2 --bfile <data> "
                    "--alt1-allele <effect_file>.snp.eff 4 3 "
                    "--make-bed --out <aligned>",
                    effect_meta.id,
                    effect_meta.A1,
                    effect_meta.A2,
                    bim_meta.A1,
                    bim_meta.A2));
        }
    }
}

auto GenotypeSelector::check_missing(const PredictObserver& observer) -> void
{
    if (missing_indices_.empty())
    {
        return;
    }

    const double missing_fraction = static_cast<double>(missing_indices_.size())
                                    / static_cast<double>(snp_effects_->size());

    if (missing_fraction > kMaxMissingFraction)
    {
        throw InvalidInputException(
            std::format(
                "Too many missing SNPs in prediction data: {}/{} ({:.1f}%)",
                missing_indices_.size(),
                snp_effects_->size(),
                missing_fraction * 100.0));
    }

    const size_t num_matched = snp_effects_->size() - missing_indices_.size();
    notify(
        observer,
        PredictSnpSelectionEvent{
            .num_matched = num_matched,
            .num_missing = missing_indices_.size(),
            .num_total = snp_effects_->size()});
}

auto GenotypeSelector::select(Eigen::MatrixXd&& raw_genotype) const
    -> Eigen::MatrixXd
{
    Eigen::MatrixXd genotype = std::move(raw_genotype);
    const Eigen::Index num_samples = genotype.rows();

    Eigen::MatrixXd result(num_samples, num_effect_snps_);

#pragma omp parallel for schedule(static) default(none) shared(result, genotype)
    for (Eigen::Index i = 0; i < num_effect_snps_; ++i)
    {
        if (column_map_[static_cast<size_t>(i)] == -1)
        {
            // Mean imputation: 2 * A1 frequency
            result.col(i).setConstant(2.0 * snp_effects_->frequencies()(i));
        }
        else
        {
            result.col(i) = genotype.col(column_map_[static_cast<size_t>(i)]);
        }
    }

    return result;
}

}  // namespace gelex
