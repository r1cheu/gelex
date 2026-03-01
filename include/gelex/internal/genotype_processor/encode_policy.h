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

#ifndef GELEX_INTERNAL_GENOTYPE_PROCESSOR_ENCODE_POLICY_H_
#define GELEX_INTERNAL_GENOTYPE_PROCESSOR_ENCODE_POLICY_H_

#include <Eigen/Core>

#include "gelex/types/genetic_effect_type.h"

namespace gelex::detail
{
template <typename T>
concept EncodePolicy
    = requires(T policy, Eigen::Ref<Eigen::VectorXd>& locus, double maf) {
          { T::encode(locus, maf) } -> std::same_as<void>;
      };

template <gelex::GeneticEffectType GT>
struct RawPolicy
{
    static auto encode(Eigen::Ref<Eigen::VectorXd> locus, double /*maf*/)
        -> void
    {
        if constexpr (GT == gelex::GeneticEffectType::Add)
        {
        }
        else
        {
            locus = locus.unaryExpr(
                [](double element) -> double
                {
                    if (element == 2.0)
                    {
                        return 0;
                    }
                    return element;
                });
        }
    };
};

template <gelex::GeneticEffectType GT>
struct OrthogonalPolicy
{
    static auto encode(Eigen::Ref<Eigen::VectorXd> locus, double maf) -> void
    {
        if constexpr (GT == gelex::GeneticEffectType::Add)
        {
        }
        else
        {
            locus = locus.unaryExpr(
                [maf](double element) -> double
                {
                    if (element == 2.0)
                    {
                        return (4 * maf) - 2;
                    }
                    if (element == 1.0)
                    {
                        return 2 * maf;
                    }
                    return element;
                });
        }
    };
};
}  // namespace gelex::detail
#endif  // GELEX_INTERNAL_GENOTYPE_PROCESSOR_ENCODE_POLICY_H_
