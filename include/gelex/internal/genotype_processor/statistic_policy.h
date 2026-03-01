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

#ifndef GELEX_INTERNAL_GENOTYPE_PROCESSOR_STATISTIC_POLICY_H_
#define GELEX_INTERNAL_GENOTYPE_PROCESSOR_STATISTIC_POLICY_H_

#include <Eigen/Core>

#include "gelex/types/genetic_effect_type.h"

namespace gelex
{
namespace detail
{
struct LocusContext
{
    Eigen::Ref<Eigen::VectorXd> locus;
    Eigen::ArrayX<bool> nan_mask;
    double valid_count;
};

template <typename T>
concept StatisticPolicy
    = requires(T policy, LocusContext& locus, LocusStatistic& statistic) {
          { T::process(locus, statistic) } -> std::same_as<void>;
      };

constexpr double MONOMORPHIC_TOL = 1e-10;

template <GeneticEffectType GT>
struct SamplePolicy
{
    static auto process(LocusContext& context, LocusStatistic& statistic)
        -> void
    {
        center(context, statistic);
        stddev(context, statistic);
    }

   private:
    static auto center(LocusContext& context, LocusStatistic& statistic) -> void
    {
        if constexpr (GT == GeneticEffectType::Add)
        {
            statistic.mean = statistic.maf * 2.0;
        }
        else
        {
            statistic.mean
                = context.nan_mask.select(0.0, context.locus.array()).sum()
                  / context.valid_count;
        }
        context.locus
            = context.nan_mask.select(statistic.mean, context.locus).array()
              - statistic.mean;
    }
    static auto stddev(LocusContext& context, LocusStatistic& statistic) -> void
    {
        statistic.stddev = std::sqrt(
            context.locus.array().square().sum() / (context.valid_count - 1.0));
        if (statistic.stddev < MONOMORPHIC_TOL)
        {
            statistic.is_monomorphic = true;
        }
    }
};

template <GeneticEffectType g_type>
struct HWEPolicy
{
    static auto process(LocusContext& context, LocusStatistic& statistic)
        -> void
    {
        center(context, statistic);
        stddev(context, statistic);
    }

   private:
    static auto center(LocusContext& context, LocusStatistic& statistic) -> void
    {
        if constexpr (g_type == GeneticEffectType::Add)
        {
            statistic.mean = statistic.maf * 2.0;
        }
        else
        {
            statistic.mean = 2 * (1 - statistic.maf) * statistic.maf;
        }
        context.locus
            = context.nan_mask.select(statistic.mean, context.locus).array()
              - statistic.mean;
    };

    static auto stddev(LocusContext& /*context*/, LocusStatistic& statistic)
        -> void
    {
        if constexpr (g_type == GeneticEffectType::Add)
        {
            statistic.stddev
                = std::sqrt(2.0 * statistic.maf * (1.0 - statistic.maf));
        }
        else
        {
            double q = 1 - statistic.maf;
            double p_sq = statistic.maf * statistic.maf;
            double q_sq = q * q;
            statistic.stddev
                = std::sqrt(2.0 * statistic.maf * q * (p_sq + q_sq));
        }
        if (statistic.stddev < MONOMORPHIC_TOL)
        {
            statistic.is_monomorphic = true;
        }
    }
};

template <GeneticEffectType g_type>
struct OrthHWEPolicy
{
    static auto process(LocusContext& context, LocusStatistic& statistic)
        -> void
    {
        center(context, statistic);
        stddev(context, statistic);
    }

   private:
    static auto center(LocusContext& context, LocusStatistic& statistic) -> void
    {
        if constexpr (g_type == GeneticEffectType::Add)
        {
            statistic.mean = statistic.maf * 2.0;
        }
        else
        {
            statistic.mean = 2.0 * statistic.maf * statistic.maf;
        }
        context.locus
            = context.nan_mask.select(statistic.mean, context.locus).array()
              - statistic.mean;
    };

    static auto stddev(LocusContext& /*context*/, LocusStatistic& statistic)
        -> void
    {
        double dominance_stddev = 2.0 * statistic.maf * (1.0 - statistic.maf);

        if constexpr (g_type == GeneticEffectType::Add)
        {
            statistic.stddev = std::sqrt(dominance_stddev);
        }
        else
        {
            statistic.stddev = dominance_stddev;
        }

        if (statistic.stddev < MONOMORPHIC_TOL)
        {
            statistic.is_monomorphic = true;
        }
    }
};

}  // namespace detail
}  // namespace gelex

#endif  // GELEX_INTERNAL_GENOTYPE_PROCESSOR_STATISTIC_POLICY_H_
