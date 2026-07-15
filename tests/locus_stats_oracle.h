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

#ifndef GELEX_TEST_LOCUS_STATS_ORACLE_H
#define GELEX_TEST_LOCUS_STATS_ORACLE_H

#include <Eigen/Core>
#include <cmath>
#include <concepts>

#include "gelex/data/encode/stats.h"

namespace gelex::test
{

// Dosage-indexed reference tabulation of LocusStats, independent of the packed
// bit-plane kernel. Serves two test roles: the oracle that count_genotypes and
// the fused encoding paths are cross-checked against, and a reliable way to
// turn a readable dosage column into the stats make_locus_encoding consumes.
template <std::floating_point T>
auto compute_locus_stats(const Eigen::Ref<const Eigen::VectorX<T>>& locus)
    -> LocusStats
{
    LocusStats stats;

    for (Eigen::Index i = 0; i < locus.size(); ++i)
    {
        const T genotype = locus[i];

        if (std::isnan(genotype))
        {
            ++stats.n_missing;
            continue;
        }

        if (genotype == T{0})
        {
            ++stats.nA2A2;
        }
        else if (genotype == T{1})
        {
            ++stats.nA1A2;
        }
        else if (genotype == T{2})
        {
            ++stats.nA1A1;
        }
    }

    return stats;
}

}  // namespace gelex::test

#endif  // GELEX_TEST_LOCUS_STATS_ORACLE_H
