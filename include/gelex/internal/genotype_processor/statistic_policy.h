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
#include "gelex/types/genotype_process_method.h"

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

template <GeneticMode GT>
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
        if constexpr (GT == GeneticMode::A)
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

template <GeneticMode GT>
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
        if constexpr (GT == GeneticMode::A)
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
        if constexpr (GT == GeneticMode::A)
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

template <GeneticMode GT>
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
        if constexpr (GT == GeneticMode::A)
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

        if constexpr (GT == GeneticMode::A)
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

template <GeneticMode GT>
struct NOIAPolicy
{
    static auto process(LocusContext& context, LocusStatistic& statistic)
        -> void
    {
        auto freqs = compute_freqs(context);
        if constexpr (GT == GeneticMode::A)
        {
            center_additive(context, freqs, statistic);
        }
        else
        {
            center_and_orthogonalize_dominance(context, freqs, statistic);
        }
        stddev(context, statistic);
    }

   private:
    struct GenotypeFreqs
    {
        double p_AA;
        double p_Aa;
        double p_aa;
    };

    static auto compute_freqs(const LocusContext& context) -> GenotypeFreqs
    {
        double n_AA
            = context.nan_mask
                  .select(0.0, (context.locus.array() == 2.0).cast<double>())
                  .sum();
        double n_Aa
            = context.nan_mask
                  .select(0.0, (context.locus.array() == 1.0).cast<double>())
                  .sum();
        double n_aa
            = context.nan_mask
                  .select(0.0, (context.locus.array() == 0.0).cast<double>())
                  .sum();
        return {
            n_AA / context.valid_count,
            n_Aa / context.valid_count,
            n_aa / context.valid_count};
    }

    static auto center_additive(
        LocusContext& context,
        const GenotypeFreqs& freqs,
        LocusStatistic& statistic) -> void
    {
        statistic.mean = (2.0 * freqs.p_AA) + freqs.p_Aa;
        context.locus
            = context.nan_mask.select(statistic.mean, context.locus).array()
              - statistic.mean;
    }

    static auto center_and_orthogonalize_dominance(
        LocusContext& context,
        const GenotypeFreqs& freqs,
        LocusStatistic& statistic) -> void
    {
        statistic.mean = 0.0;

        double diff = freqs.p_AA - freqs.p_aa;
        double denom = freqs.p_AA + freqs.p_aa - (diff * diff);
        if (denom < MONOMORPHIC_TOL)
        {
            context.locus.setZero();
            return;
        }

        double c_AA = -2.0 * freqs.p_aa * freqs.p_Aa / denom;
        double c_Aa = 4.0 * freqs.p_AA * freqs.p_aa / denom;
        double c_aa = -2.0 * freqs.p_AA * freqs.p_Aa / denom;

        context.locus = context.locus.unaryExpr(
            [c_AA, c_Aa, c_aa](double element) -> double
            {
                if (element == 2.0)
                {
                    return c_AA;
                }
                if (element == 1.0)
                {
                    return c_Aa;
                }
                if (element == 0.0)
                {
                    return c_aa;
                }
                return 0.0;
            });
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

}  // namespace detail
}  // namespace gelex

#endif  // GELEX_INTERNAL_GENOTYPE_PROCESSOR_STATISTIC_POLICY_H_
