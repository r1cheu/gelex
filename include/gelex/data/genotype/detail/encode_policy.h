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

#ifndef GELEX_DATA_GENOTYPE_DETAIL_ENCODE_POLICY_H_
#define GELEX_DATA_GENOTYPE_DETAIL_ENCODE_POLICY_H_

#include <Eigen/Core>

#include "gelex/types/genetic_effect_type.h"

namespace gelex::genotype::detail
{
template <typename T>
concept EncodePolicy
    = requires(T policy, Eigen::Ref<Eigen::VectorXd>& locus, double maf) {
          { T::encode(locus, maf) } -> std::same_as<void>;
      };

template <gelex::GeneticMode GT>
struct RawPolicy
{
    static auto encode(Eigen::Ref<Eigen::VectorXd> locus, double /*maf*/)
        -> void
    {
        if constexpr (GT == gelex::GeneticMode::A)
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

template <gelex::GeneticMode GT>
struct OrthogonalPolicy
{
    static auto encode(Eigen::Ref<Eigen::VectorXd> locus, double maf) -> void
    {
        if constexpr (GT == gelex::GeneticMode::A)
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

template <gelex::GeneticMode>
struct IdentityPolicy
{
    static auto encode(Eigen::Ref<Eigen::VectorXd> /*locus*/, double /*maf*/)
        -> void
    {
    }
};
}  // namespace gelex::genotype::detail
#endif  // GELEX_DATA_GENOTYPE_DETAIL_ENCODE_POLICY_H_
