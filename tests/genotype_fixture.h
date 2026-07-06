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

#ifndef GELEX_TESTS_GENOTYPE_FIXTURE_H_
#define GELEX_TESTS_GENOTYPE_FIXTURE_H_

#include <algorithm>
#include <cstdint>
#include <numeric>
#include <optional>
#include <utility>
#include <vector>

#include <Eigen/Core>

#include "gelex/data/genotype.h"
#include "gelex/data/snp_stats.h"

namespace gelex::test
{

class GenotypeBuilder
{
   public:
    static auto build(
        Eigen::MatrixXd data,
        std::optional<std::vector<int64_t>> valid_indices = std::nullopt,
        Eigen::VectorXd A1freq = {}) -> Genotype
    {
        OwnedStorage owned;
        owned.data = std::move(data);

        SnpStats stats;
        if (A1freq.size() == 0)
        {
            A1freq = (owned.data.colwise().mean().transpose().array() / 2.0)
                         .matrix();
        }
        stats.A1freq = std::move(A1freq);
        if (valid_indices)
        {
            stats.valid_indices = std::move(*valid_indices);
        }
        else
        {
            stats.valid_indices.resize(static_cast<size_t>(owned.data.cols()));
            std::iota(
                stats.valid_indices.begin(), stats.valid_indices.end(), 0);
        }
        std::ranges::sort(stats.valid_indices);

        owned.stats = std::move(stats);
        return Genotype{std::move(owned)};
    }
};

}  // namespace gelex::test

#endif  // GELEX_TESTS_GENOTYPE_FIXTURE_H_
