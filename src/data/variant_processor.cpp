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

#include "gelex/data/variant_processor.h"

#include <cmath>
#include <format>
#include <limits>

#include <Eigen/Dense>

#include "gelex/exception.h"

namespace gelex::detail
{

VariantStats impute_and_stats(Eigen::Ref<Eigen::VectorXd> variant)
{
    const auto n = variant.size();
    if (n < 2)
    {
        throw InvalidInputException(
            std::format("variant size {} too small for processing", n));
    }

    auto is_nan = variant.array().isNaN();
    Eigen::Index valid_count = n - is_nan.count();

    if (valid_count == 0)
    {
        variant.setZero();
        return {
            .mean = std::numeric_limits<double>::quiet_NaN(),
            .stddev = std::numeric_limits<double>::quiet_NaN(),
            .is_monomorphic = true};
    }

    double sum = is_nan.select(0.0, variant).sum();
    double mean = sum / static_cast<double>(valid_count);

    // Impute NaNs with mean
    variant = is_nan.select(mean, variant);

    const double sum_sq_diff = (variant.array() - mean).square().sum();

    double stddev = 0.0;
    if (valid_count > 1)
    {
        stddev = std::sqrt(sum_sq_diff / static_cast<double>(valid_count - 1));
    }

    return {.mean = mean, .stddev = stddev, .is_monomorphic = stddev < 1e-10};
}

VariantStats compute_and_standardize(Eigen::Ref<Eigen::VectorXd> variant)
{
    VariantStats stats = impute_and_stats(variant);

    if (!stats.is_monomorphic)
    {
        variant.array() = (variant.array() - stats.mean) / stats.stddev;
    }

    return stats;
}

GenotypeCounts count_frequencies(
    const Eigen::Ref<const Eigen::VectorXd>& variant)
{
    auto is_valid = (variant.array() == variant.array());
    double valid_count = static_cast<double>(is_valid.count());

    if (valid_count == 0)
    {
        return {0.0, 0.0, 0.0};
    }

    double count2 = static_cast<double>((variant.array() == 2.0).count());
    double count1 = static_cast<double>((variant.array() == 1.0).count());
    double count0 = valid_count - count1 - count2;

    return {
        .p_AA = count2 / valid_count,
        .p_Aa = count1 / valid_count,
        .p_aa = count0 / valid_count};
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
    return detail::impute_and_stats(variant);
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

    auto is_nan = variant.array().isNaN();
    Eigen::Index valid_count = n - is_nan.count();

    if (valid_count == 0)
    {
        VariantStats stats;
        stats.mean = std::numeric_limits<double>::quiet_NaN();
        stats.stddev = std::numeric_limits<double>::quiet_NaN();
        stats.is_monomorphic = true;
        variant.setZero();
        return stats;
    }

    double sum = is_nan.select(0.0, variant).sum();
    VariantStats stats;
    stats.mean = sum / static_cast<double>(valid_count);

    variant = is_nan.select(stats.mean, variant);

    stats.stddev = std::sqrt(stats.mean * (1.0 - 0.5 * stats.mean));
    stats.is_monomorphic = stats.stddev < 1e-10;

    if (!stats.is_monomorphic)
    {
        variant.array() = (variant.array() - stats.mean) / stats.stddev;
    }

    return stats;
}

VariantStats NOIAProcessor::process_variant(
    Eigen::Ref<Eigen::VectorXd> /*variant*/)
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

    auto is_nan = variant.array().isNaN();
    Eigen::Index valid_count = n - is_nan.count();

    if (valid_count == 0)
    {
        VariantStats stats;
        stats.mean = std::numeric_limits<double>::quiet_NaN();
        stats.stddev = std::numeric_limits<double>::quiet_NaN();
        stats.is_monomorphic = true;
        variant.setZero();
        return stats;
    }

    double sum = is_nan.select(0.0, variant).sum();
    double mean = sum / static_cast<double>(valid_count);
    const double p_freq = mean / 2.0;

    VariantStats stats;
    stats.mean = 2.0 * p_freq * p_freq;
    stats.stddev = 2.0 * p_freq * (1.0 - p_freq);
    stats.is_monomorphic = stats.stddev < 1e-10;

    const double one_alt_encode = 2.0 * p_freq;
    const double two_alt_encode = (4.0 * p_freq) - 2.0;

    variant = variant.unaryExpr(
        [&](double x) -> double
        {
            if (std::isnan(x))
                return mean;
            if (x == 1.0)
                return one_alt_encode;
            if (x == 2.0)
                return two_alt_encode;
            return x;
        });

    if (!stats.is_monomorphic)
    {
        variant.array() = (variant.array() - stats.mean) / stats.stddev;
    }

    return stats;
}

VariantStats DominantNOIAProcessor::process_variant(
    Eigen::Ref<Eigen::VectorXd> /*variant*/)
{
    throw InvalidInputException(
        "Dominant NOIA processing is not implemented yet.");
}

}  // namespace gelex
