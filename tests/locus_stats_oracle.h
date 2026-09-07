// Copyright 2026 RuLei Chen
// SPDX-License-Identifier: Apache-2.0

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
