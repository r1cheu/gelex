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

#include "gelex/predict/standardize.h"

#include <cmath>

#include <Eigen/Core>

#include "gelex/data/genotype_method.h"
#include "gelex/predict/types.h"

namespace gelex
{

auto standardize_genotypes(GenotypeData& geno, const SnpStatsData& snpstats)
    -> void
{
    const auto n_snps = geno.add.cols();
    const bool add_scale = !is_center(snpstats.add.method);

    for (Eigen::Index j = 0; j < n_snps; ++j)
    {
        auto col = geno.add.col(j);
        col.array() -= snpstats.add.mean(j);
        if (add_scale)
        {
            col.array() /= std::sqrt(snpstats.add.var(j));
        }
        col = col.array().isNaN().select(0.0, col.array());
    }

    if (geno.dom.has_value())
    {
        const bool use_orthogonal = is_orthogonal(snpstats.add.method);
        const bool dom_scale = !is_center(snpstats.dom.method);

        for (Eigen::Index j = 0; j < n_snps; ++j)
        {
            auto col = geno.dom->col(j);
            const double maf = snpstats.add.mean(j) / 2.0;

            if (use_orthogonal)
            {
                col = col.unaryExpr(
                    [maf](double genotype) -> double
                    {
                        if (genotype == 2.0)
                        {
                            return (4.0 * maf) - 2.0;
                        }
                        if (genotype == 1.0)
                        {
                            return 2.0 * maf;
                        }
                        return genotype;
                    });
            }
            else
            {
                col = col.unaryExpr(
                    [](double genotype) -> double
                    {
                        if (genotype == 2.0)
                        {
                            return 0.0;
                        }
                        return genotype;
                    });
            }

            col.array() -= snpstats.dom.mean(j);
            if (dom_scale)
            {
                col.array() /= std::sqrt(snpstats.dom.var(j));
            }
            col = col.array().isNaN().select(0.0, col.array());
        }
    }
}

}  // namespace gelex
