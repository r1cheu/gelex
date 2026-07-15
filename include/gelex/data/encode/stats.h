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
