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

#ifndef GELEX_DATA_GENOTYPE_PROCESSOR_H_
#define GELEX_DATA_GENOTYPE_PROCESSOR_H_

#include <algorithm>
#include <cmath>
#include <concepts>

#include <omp.h>
#include <Eigen/Dense>

namespace gelex
{

constexpr double GENOTYPE_EPSILON = 1e-10;

struct VariantStats
{
    double mean{0.0};
    double stddev{0.0};
    bool is_monomorphic{false};
};

template <typename T>
concept GenotypeProcessor
    = requires(T processor, Eigen::Ref<Eigen::VectorXd> variant) {
          { T::process_variant(variant) } -> std::same_as<VariantStats>;
          { T::dom } -> std::convertible_to<bool>;
      };

namespace detail
{

using NanMask = Eigen::Array<bool, Eigen::Dynamic, 1>;

struct VariantContext
{
    NanMask is_nan;
    Eigen::Index valid_count{0};
    double additive_sample_mean{0.0};
    double p{0.0};
};

inline auto compute_mean(Eigen::Ref<const Eigen::VectorXd> variant) -> double
{
    auto is_nan = variant.array().isNaN();
    Eigen::Index valid_count = variant.size() - is_nan.count();
    if (valid_count <= 0)
    {
        return 0.0;
    }
    double sum = is_nan.select(0.0, variant).sum();
    return sum / static_cast<double>(valid_count);
}

inline auto compute_sample_stddev(
    Eigen::Ref<const Eigen::VectorXd> centered,
    Eigen::Index valid_count) -> double
{
    if (valid_count <= 1)
    {
        return 0.0;
    }
    double sum_sq_diff = centered.array().square().sum();
    return std::sqrt(sum_sq_diff / (static_cast<double>(valid_count) - 1.0));
}

inline auto estimate_p_from_additive(double additive_mean) -> double
{
    return std::clamp(additive_mean / 2.0, 0.0, 1.0);
}

inline auto additive_hwe_mean(double p) -> double
{
    return 2.0 * p;
}

inline auto additive_hwe_stddev(double p) -> double
{
    double q = 1.0 - p;
    return std::sqrt(std::max(0.0, 2.0 * p * q));
}

inline auto build_context(Eigen::Ref<const Eigen::VectorXd> variant)
    -> VariantContext
{
    VariantContext context;
    context.is_nan = variant.array().isNaN();
    context.valid_count = variant.size() - context.is_nan.count();

    if (context.valid_count <= 0)
    {
        return context;
    }

    context.additive_sample_mean = context.is_nan.select(0.0, variant).sum()
                                   / static_cast<double>(context.valid_count);
    context.p = estimate_p_from_additive(context.additive_sample_mean);
    return context;
}

}  // namespace detail

struct DomBinary
{
    static auto encode(double genotype, double /*p*/) -> double
    {
        if (genotype == 2.0)
        {
            return 0.0;
        }
        return genotype;
    }
};

struct DomOrthogonal
{
    static auto encode(double genotype, double p) -> double
    {
        if (genotype == 2.0)
        {
            return (4.0 * p) - 2.0;
        }
        if (genotype == 1.0)
        {
            return 2.0 * p;
        }
        return 0.0;
    }
};

struct StatsSample
{
    static auto additive_mean(double additive_sample_mean, double /*p*/)
        -> double
    {
        return additive_sample_mean;
    }

    static auto additive_stddev(
        Eigen::Ref<const Eigen::VectorXd> centered,
        Eigen::Index valid_count,
        double /*p*/) -> double
    {
        return detail::compute_sample_stddev(centered, valid_count);
    }

    static auto dominant_mean(double dominant_sample_mean, double /*p*/)
        -> double
    {
        return dominant_sample_mean;
    }

    static auto dominant_stddev(
        Eigen::Ref<const Eigen::VectorXd> centered,
        Eigen::Index valid_count,
        double /*p*/) -> double
    {
        return detail::compute_sample_stddev(centered, valid_count);
    }
};

struct StatsHWE
{
    static auto additive_mean(double /*additive_sample_mean*/, double p)
        -> double
    {
        return detail::additive_hwe_mean(p);
    }

    static auto additive_stddev(
        Eigen::Ref<const Eigen::VectorXd> /*centered*/,
        Eigen::Index /*valid_count*/,
        double p) -> double
    {
        return detail::additive_hwe_stddev(p);
    }

    static auto dominant_mean(double /*dominant_sample_mean*/, double p)
        -> double
    {
        double q = 1.0 - p;
        return 2.0 * p * q;
    }

    static auto dominant_stddev(
        Eigen::Ref<const Eigen::VectorXd> /*centered*/,
        Eigen::Index /*valid_count*/,
        double p) -> double
    {
        double q = 1.0 - p;
        double p_sq = p * p;
        double q_sq = q * q;
        return std::sqrt(std::max(0.0, 2.0 * p * q * (p_sq + q_sq)));
    }
};

struct StatsOrthHWE
{
    static auto additive_mean(double /*additive_sample_mean*/, double p)
        -> double
    {
        return detail::additive_hwe_mean(p);
    }

    static auto additive_stddev(
        Eigen::Ref<const Eigen::VectorXd> /*centered*/,
        Eigen::Index /*valid_count*/,
        double p) -> double
    {
        return detail::additive_hwe_stddev(p);
    }

    static auto dominant_mean(double /*dominant_sample_mean*/, double p)
        -> double
    {
        return 2.0 * p * p;
    }

    static auto dominant_stddev(
        Eigen::Ref<const Eigen::VectorXd> /*centered*/,
        Eigen::Index /*valid_count*/,
        double p) -> double
    {
        double q = 1.0 - p;
        return 2.0 * p * q;
    }
};

struct ScaleNone
{
    static auto apply(
        Eigen::Ref<Eigen::VectorXd> /*variant*/,
        double /*stddev*/) -> void
    {
    }
};

struct ScaleStandardize
{
    static auto apply(Eigen::Ref<Eigen::VectorXd> variant, double stddev)
        -> void
    {
        variant.array() /= stddev;
    }
};

template <
    typename DominantEncodingPolicy,
    typename StatsPolicy,
    typename ScalingPolicy>
struct GenotypeStrategy
{
   private:
    static auto set_missing_and_center(
        Eigen::Ref<Eigen::VectorXd> variant,
        const detail::NanMask& is_nan,
        double mean) -> void
    {
        variant = is_nan.select(mean, variant);
        variant.array() -= mean;
    }

    static auto finalize_variant(
        Eigen::Ref<Eigen::VectorXd> variant,
        double mean,
        double stddev) -> VariantStats
    {
        VariantStats stats{mean, stddev, false};
        if (stddev < GENOTYPE_EPSILON)
        {
            stats.is_monomorphic = true;
            return stats;
        }

        ScalingPolicy::apply(variant, stddev);
        return stats;
    }

   public:
    static auto process_additive(Eigen::Ref<Eigen::VectorXd> variant)
        -> VariantStats
    {
        auto context = detail::build_context(variant);
        if (context.valid_count <= 0)
        {
            variant.setZero();
            return {0.0, 0.0, true};
        }

        double mean = StatsPolicy::additive_mean(
            context.additive_sample_mean, context.p);
        set_missing_and_center(variant, context.is_nan, mean);

        double stddev = StatsPolicy::additive_stddev(
            variant, context.valid_count, context.p);

        return finalize_variant(variant, mean, stddev);
    }

    static auto process_dominant(Eigen::Ref<Eigen::VectorXd> variant)
        -> VariantStats
    {
        auto context = detail::build_context(variant);
        if (context.valid_count <= 0)
        {
            variant.setZero();
            return {0.0, 0.0, true};
        }

        variant = variant.unaryExpr(
            [p = context.p](double genotype) -> double
            {
                if (std::isnan(genotype))
                {
                    return genotype;
                }
                return DominantEncodingPolicy::encode(genotype, p);
            });

        double dominant_sample_mean = detail::compute_mean(variant);
        double mean
            = StatsPolicy::dominant_mean(dominant_sample_mean, context.p);

        set_missing_and_center(variant, context.is_nan, mean);

        double stddev = StatsPolicy::dominant_stddev(
            variant, context.valid_count, context.p);

        return finalize_variant(variant, mean, stddev);
    }
};

using CenterMethod = GenotypeStrategy<DomBinary, StatsSample, ScaleNone>;
using StandardizeMethod
    = GenotypeStrategy<DomBinary, StatsSample, ScaleStandardize>;
using OrthCenterMethod
    = GenotypeStrategy<DomOrthogonal, StatsSample, ScaleNone>;
using OrthStandardizeMethod
    = GenotypeStrategy<DomOrthogonal, StatsSample, ScaleStandardize>;

using CenterHWEMethod = GenotypeStrategy<DomBinary, StatsHWE, ScaleNone>;
using StandardizeHWEMethod
    = GenotypeStrategy<DomBinary, StatsHWE, ScaleStandardize>;
using OrthCenterHWEMethod
    = GenotypeStrategy<DomOrthogonal, StatsOrthHWE, ScaleNone>;
using OrthStandardizeHWEMethod
    = GenotypeStrategy<DomOrthogonal, StatsOrthHWE, ScaleStandardize>;

// =========================================================================
// Template adapters: adapt Method traits to GenotypeProcessor concept
// =========================================================================

template <typename Method>
struct AdditiveProcessor
{
    static constexpr bool dom = false;
    static auto process_variant(Eigen::Ref<Eigen::VectorXd> variant)
        -> VariantStats
    {
        return Method::process_additive(variant);
    }
};

template <typename Method>
struct DominantProcessor
{
    static constexpr bool dom = true;
    static auto process_variant(Eigen::Ref<Eigen::VectorXd> variant)
        -> VariantStats
    {
        return Method::process_dominant(variant);
    }
};

// =========================================================================
// Batch processing template
// =========================================================================

template <GenotypeProcessor P>
auto process_matrix(
    Eigen::Ref<Eigen::MatrixXd> genotype,
    Eigen::VectorXd* freqs = nullptr) -> void
{
#pragma omp parallel for default(none) shared(genotype, freqs)
    for (Eigen::Index i = 0; i < genotype.cols(); ++i)
    {
        auto col = genotype.col(i);
        auto stats = P::process_variant(col);
        if (freqs != nullptr)
        {
            (*freqs)(i) = stats.mean / 2;
        }
    }
}

// =========================================================================
// GRM method bundles
// =========================================================================

namespace grm
{

struct Standardized
{
    using Additive = AdditiveProcessor<StandardizeMethod>;
    using Dominant = DominantProcessor<StandardizeMethod>;
};

struct OrthStandardized
{
    using Additive = AdditiveProcessor<OrthStandardizeMethod>;
    using Dominant = DominantProcessor<OrthStandardizeMethod>;
};

struct Centered
{
    using Additive = AdditiveProcessor<CenterMethod>;
    using Dominant = DominantProcessor<CenterMethod>;
};

struct OrthCentered
{
    using Additive = AdditiveProcessor<OrthCenterMethod>;
    using Dominant = DominantProcessor<OrthCenterMethod>;
};

struct StandardizedHWE
{
    using Additive = AdditiveProcessor<StandardizeHWEMethod>;
    using Dominant = DominantProcessor<StandardizeHWEMethod>;
};

struct OrthStandardizedHWE
{
    using Additive = AdditiveProcessor<OrthStandardizeHWEMethod>;
    using Dominant = DominantProcessor<OrthStandardizeHWEMethod>;
};

struct CenteredHWE
{
    using Additive = AdditiveProcessor<CenterHWEMethod>;
    using Dominant = DominantProcessor<CenterHWEMethod>;
};

struct OrthCenteredHWE
{
    using Additive = AdditiveProcessor<OrthCenterHWEMethod>;
    using Dominant = DominantProcessor<OrthCenterHWEMethod>;
};

}  // namespace grm

}  // namespace gelex

#endif  // GELEX_DATA_GENOTYPE_PROCESSOR_H_
