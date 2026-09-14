// Copyright 2026 RuLei Chen
// SPDX-License-Identifier: Apache-2.0

#include "gelex/bayes/genotype/compact_genotype.h"

#include <Eigen/Core>
#include <cstddef>
#include <cstdint>
#include <fmt/format.h>
#include <functional>
#include <span>
#include <utility>
#include <vector>

#include "gelex/data/bed.h"
#include "gelex/data/encode/encoder.h"
#include "gelex/data/encode/stats.h"
#include "gelex/exception.h"
#include "gelex/infra/notify.h"

namespace gelex::bayes
{

CompactGenotype::CompactGenotype(
    raw_matrix_type raw_codes,
    std::vector<gelex::LocusStats> locus_stats,
    Eigen::VectorXd a1_frequency)
    : raw_codes_{std::move(raw_codes)},
      locus_stats_{std::move(locus_stats)},
      a1_frequency_{std::move(a1_frequency)}
{
    const auto markers = static_cast<std::size_t>(raw_codes_.cols());
    if (locus_stats_.size() != markers
        || static_cast<std::size_t>(a1_frequency_.size()) != markers)
    {
        throw GelexException{fmt::format(
            "CompactGenotype: raw codes have {} markers but locus stats have "
            "{} and A1 frequencies have {}",
            markers,
            locus_stats_.size(),
            a1_frequency_.size())};
    }
}

auto make_compact_genotype(
    const gelex::Bed& bed,
    const std::function<void(std::size_t)>& observer) -> CompactGenotype
{
    const auto rows = bed.num_samples();
    const auto cols = bed.num_snps();
    CompactGenotype::raw_matrix_type raw_codes(rows, cols);
    std::vector<gelex::LocusStats> locus_stats(static_cast<std::size_t>(cols));
    Eigen::VectorXd a1_frequency(cols);

    const gelex::LocusEncoder encoder{bed};
    for (Eigen::Index marker = 0; marker < cols; ++marker)
    {
        encoder.decode_raw_codes(
            marker,
            std::span<std::uint8_t>{
                raw_codes.col(marker).data(), static_cast<std::size_t>(rows)});
        const auto stats = encoder.count(marker);
        locus_stats[static_cast<std::size_t>(marker)] = stats;
        a1_frequency[marker] = stats.has_nonmissing() ? stats.A1freq() : 0.0;
        notify(observer, static_cast<std::size_t>(marker + 1));
    }
    return CompactGenotype{
        std::move(raw_codes), std::move(locus_stats), std::move(a1_frequency)};
}

}  // namespace gelex::bayes
