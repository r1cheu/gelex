// Copyright 2026 RuLei Chen
// SPDX-License-Identifier: Apache-2.0

#ifndef GELEX_DATA_ENCODE_STATS_H_
#define GELEX_DATA_ENCODE_STATS_H_

#include <Eigen/Core>

namespace gelex
{

struct LocusStats
{
    Eigen::Index nA2A2{0};
    Eigen::Index nA1A2{0};
    Eigen::Index nA1A1{0};
    Eigen::Index n_missing{0};

    [[nodiscard]] auto n_nonmissing() const noexcept -> Eigen::Index
    {
        return nA2A2 + nA1A2 + nA1A1;
    }

    [[nodiscard]] auto has_nonmissing() const noexcept -> bool
    {
        return n_nonmissing() > 0;
    }

    [[nodiscard]] auto pA2A2() const -> double
    {
        return static_cast<double>(nA2A2) / static_cast<double>(n_nonmissing());
    }

    [[nodiscard]] auto pA1A2() const -> double
    {
        return static_cast<double>(nA1A2) / static_cast<double>(n_nonmissing());
    }

    [[nodiscard]] auto pA1A1() const -> double
    {
        return static_cast<double>(nA1A1) / static_cast<double>(n_nonmissing());
    }

    [[nodiscard]] auto A1freq() const -> double
    {
        return pA1A1() + (0.5 * pA1A2());
    }
};

}  // namespace gelex

#endif  // GELEX_DATA_ENCODE_STATS_H_
