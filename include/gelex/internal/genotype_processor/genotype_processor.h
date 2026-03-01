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

#ifndef GELEX_INTERNAL_GENOTYPE_PROCESSOR_GENOTYPE_PROCESSOR_H_
#define GELEX_INTERNAL_GENOTYPE_PROCESSOR_GENOTYPE_PROCESSOR_H_

#include <Eigen/Dense>

#include "gelex/internal/genotype_processor/encode_policy.h"
#include "gelex/internal/genotype_processor/statistic_policy.h"
#include "gelex/types/genetic_effect_type.h"

namespace gelex::detail
{

template <
    detail::EncodePolicy Encode,
    detail::StatisticPolicy Stats,
    bool Scale>
struct GenotypeProcessorStrategy
{
    static auto process(Eigen::Ref<Eigen::VectorXd> locus) -> LocusStatistic
    {
        auto nan_mask = locus.array().isNaN();
        auto valid_count = static_cast<double>(locus.size() - nan_mask.count());
        if (valid_count <= 0)
        {
            return LocusStatistic{};
        }
        LocusStatistic stats{};
        stats.maf = nan_mask.select(0.0, locus).sum() / (2.0 * valid_count);
        LocusContext context{locus, nan_mask, valid_count};
        Encode::encode(locus, stats.maf);
        Stats::process(context, stats);

        if (stats.is_monomorphic)
        {
            return stats;
        }
        if constexpr (Scale)
        {
            locus.array() /= stats.stddev;
        }
        return stats;
    }
};

}  // namespace gelex::detail

#endif  // GELEX_INTERNAL_GENOTYPE_PROCESSOR_GENOTYPE_PROCESSOR_H_
