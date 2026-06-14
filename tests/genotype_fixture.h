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
#include <utility>
#include <vector>

#include <Eigen/Core>

#include "gelex/data/genotype/genotype.h"

namespace gelex::test
{

class GenotypeBuilder
{
   public:
    static auto build(
        Eigen::MatrixXd data,
        Eigen::VectorXd mean,
        Eigen::VectorXd var,
        std::vector<int64_t> mono_indices = {},
        Eigen::VectorXd A1freq = {}) -> genotype::Genotype
    {
        genotype::OwnedStorage owned;
        owned.data = std::move(data);
        if (A1freq.size() == 0)
        {
            A1freq = (mean.array() / 2.0).matrix();
        }
        owned.mean = std::move(mean);
        owned.var = std::move(var);
        owned.A1freq = std::move(A1freq);
        owned.mono_indices = std::move(mono_indices);
        std::ranges::sort(owned.mono_indices);
        return genotype::Genotype{std::move(owned)};
    }
};

}  // namespace gelex::test

#endif  // GELEX_TESTS_GENOTYPE_FIXTURE_H_
