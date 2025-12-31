#include "gelex/data/variant_processor.h"

#include <cmath>
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
        throw InvalidInputException(
            std::format("variant size {} too small for processing", n));
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
    auto count2 = static_cast<double>((variant.array() == 2.0).count());
    auto count1 = static_cast<double>((variant.array() == 1.0).count());
    auto n = static_cast<double>(variant.size());

    double count0 = n - count1 - count2;

    return {.p_AA = count2 / n, .p_Aa = count1 / n, .p_aa = count0 / n};
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
        throw InvalidInputException(
            std::format("variant size {} too small for processing", n));
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
        throw InvalidInputException(
            std::format("variant size {} too small for processing", n));
    }

    VariantStats stats;
    stats.mean = variant.mean();
    stats.stddev
        = std::sqrt(stats.mean * (1.0 - 0.5 * stats.mean));  // sqrt(2pq)
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
    throw InvalidInputException(
        "Dominant NOIA processing is not implemented yet.");
}

VariantStats DominantStandardizingProcessor::process_variant(
    Eigen::Ref<Eigen::VectorXd> variant)
{
    const auto n = variant.size();
    if (n < 2)
    {
        throw InvalidInputException(
            std::format("variant size {} too small for processing", n));
    }
    variant = (variant.array() == 2.0).select(0.0, variant);

    return detail::compute_and_standardize(variant);
}

VariantStats DominantRawProcessor::process_variant(
    Eigen::Ref<Eigen::VectorXd> variant)
{
    const auto n = variant.size();
    if (n < 2)
    {
        throw InvalidInputException(
            std::format("variant size {} too small for processing", n));
    }
    variant = (variant.array() == 2.0).select(0.0, variant);
    return RawProcessor::process_variant(variant);
}

VariantStats DominantOrthogonalHWEProcessor::process_variant(
    Eigen::Ref<Eigen::VectorXd> variant)
{
    const auto n = variant.size();
    if (n < 2)
    {
        throw InvalidInputException(
            std::format("variant size {} too small for processing", n));
    }
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

VariantStats DominantNOIAProcessor::process_variant(
    Eigen::Ref<Eigen::VectorXd> variant)
{
    throw InvalidInputException(
        "Dominant NOIA processing is not implemented yet.");
};
}  // namespace gelex
