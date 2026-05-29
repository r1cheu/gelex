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

#include "standardize.h"

#include <Eigen/Core>

#include "gelex/data/genotype/detail/encode_policy.h"
#include "gelex/data/genotype/process_method.h"
#include "gelex/predict/types.h"
#include "gelex/types/genetic_effect_type.h"

namespace gelex::predict::detail
{

auto standardize_genotypes(GenotypeData& geno, const SbinData& sbin) -> void
{
    const auto n_snps = geno.add.cols();
    const bool has_stddev = sbin.add.stddev.has_value();

    for (Eigen::Index j = 0; j < n_snps; ++j)
    {
        auto col = geno.add.col(j);
        col.array() -= sbin.add.mean(j);
        if (has_stddev)
        {
            col.array() /= (*sbin.add.stddev)(j);
        }
        col = col.array().isNaN().select(0.0, col.array());
    }

    if (geno.dom.has_value())
    {
        const bool use_orthogonal = sbin.add.method.is_orthogonal();
        const bool dom_has_stddev = sbin.dom.stddev.has_value();

        for (Eigen::Index j = 0; j < n_snps; ++j)
        {
            auto col = geno.dom->col(j);
            const double maf = sbin.add.mean(j) / 2.0;

            if (use_orthogonal)
            {
                gelex::genotype::detail::OrthogonalPolicy<
                    GeneticMode::D>::encode(col, maf);
            }
            else
            {
                gelex::genotype::detail::RawPolicy<GeneticMode::D>::encode(
                    col, maf);
            }

            col.array() -= sbin.dom.mean(j);
            if (dom_has_stddev)
            {
                col.array() /= (*sbin.dom.stddev)(j);
            }
            col = col.array().isNaN().select(0.0, col.array());
        }
    }
}

}  // namespace gelex::predict::detail
