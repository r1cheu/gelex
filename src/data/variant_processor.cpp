#include "gelex/data/variant_processor.h"

#include <cmath>
#include <cstddef>
#include <format>
#include <limits>

#include <Eigen/Dense>

#include "gelex/exception.h"

namespace gelex::detail
{

VariantStats compute_and_standardize(Eigen::Ref<Eigen::VectorXd> variant)
{
    const auto n = variant.size();
    if (n < 2)
    {
        throw InvalidDataException(
            std::format("Variant size too small for processing: {}", n));
    }

    VariantStats stats;
    stats.mean = variant.mean();

    const double sum_sq_diff = (variant.array() - stats.mean).square().sum();
    stats.stddev = std::sqrt(sum_sq_diff / static_cast<double>(n - 1));

    stats.is_monomorphic
        = stats.stddev < std::numeric_limits<double>::epsilon();

    if (!stats.is_monomorphic)
    {
        variant.array() = (variant.array() - stats.mean) / stats.stddev;
    }

    return stats;
}

GenotypeCounts count_frequencies(
    const Eigen::Ref<const Eigen::VectorXd>& variant)
{
    size_t count0 = 0;
    size_t count1 = 0;
    size_t count2 = 0;
    const double* data = variant.data();
    const auto size = variant.size();
    const auto stride = variant.innerStride();

    for (Eigen::Index i = 0; i < size; ++i)
    {
        double val = data[i * stride];
        if (val == 1.0)
        {
            count1++;
        }
        else if (val == 2.0)
        {
            count2++;
        }
        else
        {
            count0++;
        }
    }

    auto n = static_cast<double>(size);
    return {
        .p_AA = static_cast<double>(count2) / n,
        .p_Aa = static_cast<double>(count1) / n,
        .p_aa = static_cast<double>(count0) / n};
}

}  // namespace gelex::detail

namespace gelex
{

VariantStats StandardizingProcessor::process_variant(
    Eigen::Ref<Eigen::VectorXd> variant)
{
    return detail::compute_and_standardize(variant);
}

VariantStats RawProcessor::process_variant(Eigen::Ref<Eigen::VectorXd> variant)
{
    const auto n = variant.size();
    if (n < 2)
    {
        throw InvalidDataException("Variant size too small");
    }

    VariantStats stats;
    stats.mean = variant.mean();
    double sum_sq_diff = (variant.array() - stats.mean).square().sum();
    stats.stddev = std::sqrt(sum_sq_diff / static_cast<double>(n - 1));
    stats.is_monomorphic
        = stats.stddev < std::numeric_limits<double>::epsilon();

    return stats;
}

VariantStats HardWenbergProcessor::process_variant(
    Eigen::Ref<Eigen::VectorXd> variant)
{
    const auto n = variant.size();
    if (n < 2)
    {
        throw InvalidDataException("Variant size too small");
    }

    VariantStats stats;
    stats.mean = variant.mean();
    // HWE 假设下的标准差公式
    stats.stddev = std::sqrt(stats.mean * (1.0 - 0.5 * stats.mean));
    stats.is_monomorphic
        = stats.stddev < std::numeric_limits<double>::epsilon();

    if (!stats.is_monomorphic)
    {
        variant.array() = (variant.array() - stats.mean) / stats.stddev;
    }
    return stats;
}

VariantStats NOIAProcessor::process_variant(Eigen::Ref<Eigen::VectorXd> variant)
{
    auto counts = detail::count_frequencies(variant);

    const double AA = -(-counts.p_Aa - (2 * counts.p_aa));
    const double Aa = -(1 - counts.p_Aa - (2 * counts.p_aa));
    const double aa = -(2 - counts.p_Aa - (2 * counts.p_aa));

    variant = variant.unaryExpr(
        [&](double x) -> double
        {
            if (x == 0.0)
            {
                return AA;
            }
            if (x == 1.0)
            {
                return Aa;
            }
            return aa;
        });

    return detail::compute_and_standardize(variant);
}

VariantStats DominantStandardizingProcessor::process_variant(
    Eigen::Ref<Eigen::VectorXd> variant)
{
    variant = (variant.array() == 2.0).select(0.0, variant);

    return detail::compute_and_standardize(variant);
}

VariantStats DominantRawProcessor::process_variant(
    Eigen::Ref<Eigen::VectorXd> variant)
{
    variant = (variant.array() == 2.0).select(0.0, variant);
    return RawProcessor::process_variant(variant);
}

VariantStats DominantOrthogonalHWEProcessor::process_variant(
    Eigen::Ref<Eigen::VectorXd> variant)
{
    const double p_freq = variant.mean() / 2.0;

    VariantStats stats;
    stats.mean = 2.0 * p_freq * p_freq;
    stats.stddev = 2.0 * p_freq * (1.0 - p_freq);
    stats.is_monomorphic
        = stats.stddev < std::numeric_limits<double>::epsilon();

    const double one_alt_encode = 2.0 * p_freq;
    const double two_alt_encode = (4.0 * p_freq) - 2.0;

    variant = variant.unaryExpr(
        [&](double x) -> double
        {
            if (x == 1.0)
            {
                return one_alt_encode;
            }
            if (x == 2.0)
            {
                return two_alt_encode;
            }
            return x;
        });

    if (!stats.is_monomorphic)
    {
        variant.array() = (variant.array() - stats.mean) / stats.stddev;
    }
    return stats;
}

VariantStats DominantNOIAHWEProcessor::process_variant(
    Eigen::Ref<Eigen::VectorXd> variant)
{
    auto counts = detail::count_frequencies(variant);

    const double denom
        = counts.p_AA + counts.p_aa - std::pow(counts.p_AA - counts.p_aa, 2);

    if (std::abs(denom) < std::numeric_limits<double>::epsilon())
    {
        throw InvalidOperationException(
            "DominantNOIAHWE calculation resulted in zero denominator");
    }

    const double AA = -(2 * counts.p_Aa * counts.p_aa) / denom;
    const double Aa = (4 * counts.p_AA * counts.p_aa) / denom;
    const double aa = -(2 * counts.p_AA * counts.p_Aa) / denom;

    variant = variant.unaryExpr(
        [&](double x) -> double
        {
            if (x == 0.0)
            {
                return AA;
            }
            if (x == 1.0)
            {
                return Aa;
            }
            return aa;
        });
    return detail::compute_and_standardize(variant);
}

}  // namespace gelex
