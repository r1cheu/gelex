// Copyright 2026 RuLei Chen
// SPDX-License-Identifier: Apache-2.0

#include "gelex/data/encode/matrix.h"

#include <Eigen/Core>
#include <cmath>
#include <fmt/format.h>

#include "gelex/data/encode/detail/encoding.h"
#include "gelex/data/encode/spec.h"
#include "gelex/data/encode/stats.h"
#include "gelex/data/encode/types.h"
#include "gelex/exception.h"

namespace gelex
{

namespace
{

auto count_genotypes(
    const Eigen::Ref<const Eigen::VectorXd>& genotypes,
    Eigen::Index marker) -> LocusStats
{
    LocusStats stats;
    for (Eigen::Index sample = 0; sample < genotypes.size(); ++sample)
    {
        const double genotype = genotypes(sample);
        if (std::isnan(genotype))
        {
            ++stats.n_missing;
        }
        else if (genotype == 0.0)
        {
            ++stats.nA2A2;
        }
        else if (genotype == 1.0)
        {
            ++stats.nA1A2;
        }
        else if (genotype == 2.0)
        {
            ++stats.nA1A1;
        }
        else
        {
            throw GelexException(
                fmt::format(
                    "genotype at sample {} and marker {} must be 0, 1, 2, or "
                    "NaN",
                    sample,
                    marker));
        }
    }
    return stats;
}

auto apply_encoding(
    Eigen::Ref<Eigen::VectorXd> genotypes,
    const LocusEncoding& encoding) -> void
{
    if (!encoding.valid)
    {
        genotypes.setZero();
        return;
    }

    for (Eigen::Index sample = 0; sample < genotypes.size(); ++sample)
    {
        const double genotype = genotypes(sample);
        if (std::isnan(genotype))
        {
            genotypes(sample) = encoding.lut[1];
        }
        else if (genotype == 0.0)
        {
            genotypes(sample) = encoding.lut[3];
        }
        else if (genotype == 1.0)
        {
            genotypes(sample) = encoding.lut[2];
        }
        else
        {
            genotypes(sample) = encoding.lut[0];
        }
    }
}

}  // namespace

auto encode_inplace(
    Eigen::Ref<Eigen::MatrixXd> genotypes,
    GeneticMode mode,
    GenotypeMethod method) -> void
{
    const EncodingSpec spec{encoding_spec_from_method(mode, method)};
    for (Eigen::Index marker = 0; marker < genotypes.cols(); ++marker)
    {
        auto column = genotypes.col(marker);
        const LocusStats stats{count_genotypes(column, marker)};
        const LocusEncoding encoding{
            detail::make_locus_encoding(marker, stats, spec)};
        apply_encoding(column, encoding);
    }
}

}  // namespace gelex
