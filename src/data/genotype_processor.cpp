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

#include "gelex/data/genotype_processor.h"

#include <Eigen/Dense>

namespace gelex
{

// =========================================================================
// CenterMethod
// =========================================================================

auto CenterMethod::process_additive(Eigen::Ref<Eigen::VectorXd> variant)
    -> VariantStats
{
    auto is_nan = variant.array().isNaN();
    Eigen::Index valid_count = variant.size() - is_nan.count();

    double sum = is_nan.select(0.0, variant).sum();
    double mean = sum / static_cast<double>(valid_count);
    variant = is_nan.select(mean, variant);

    double sum_sq_diff = (variant.array() - mean).square().sum();
    double stddev
        = std::sqrt(sum_sq_diff / (static_cast<double>(valid_count) - 1));

    variant.array() -= mean;
    VariantStats stats{mean, stddev, false};
    if (stddev < GENOTYPE_EPSILON)
    {
        stats.is_monomorphic = true;
        return stats;
    }
    return stats;
}

auto CenterMethod::process_dominant(Eigen::Ref<Eigen::VectorXd> variant)
    -> VariantStats
{
    variant = (variant.array() == 2.0).select(0.0, variant);
    return process_additive(variant);
}

// =========================================================================
// StandardizeMethod
// =========================================================================

auto StandardizeMethod::process_additive(Eigen::Ref<Eigen::VectorXd> variant)
    -> VariantStats
{
    auto stats = CenterMethod::process_additive(variant);
    if (stats.is_monomorphic)
    {
        return stats;
    }
    variant.array() /= stats.stddev;
    return stats;
}

auto StandardizeMethod::process_dominant(Eigen::Ref<Eigen::VectorXd> variant)
    -> VariantStats
{
    variant = (variant.array() == 2.0).select(0.0, variant);
    return process_additive(variant);
}

// =========================================================================
// OrthCenterMethod
// =========================================================================

auto OrthCenterMethod::process_additive(Eigen::Ref<Eigen::VectorXd> variant)
    -> VariantStats
{
    return CenterMethod::process_additive(variant);
}

auto OrthCenterMethod::process_dominant(Eigen::Ref<Eigen::VectorXd> variant)
    -> VariantStats
{
    auto is_nan = variant.array().isNaN();
    Eigen::Index valid_count = variant.size() - is_nan.count();

    double sum = is_nan.select(0.0, variant).sum();
    double mean = sum / static_cast<double>(valid_count);

    double maf = mean / 2;
    double one_alt_encode = 2.0 * maf;
    double two_alt_encode = (4.0 * maf) - 2.0;

    variant = variant.unaryExpr(
        [one_alt_encode, two_alt_encode](double x) -> double
        {
            if (std::isnan(x))
            {
                return x;
            }
            if (x == 2.0)
            {
                return two_alt_encode;
            }
            if (x == 1.0)
            {
                return one_alt_encode;
            }
            return 0.0;
        });

    return CenterMethod::process_additive(variant);
}

// =========================================================================
// OrthStandardizeMethod
// =========================================================================

auto OrthStandardizeMethod::process_additive(
    Eigen::Ref<Eigen::VectorXd> variant) -> VariantStats
{
    auto stats = CenterMethod::process_additive(variant);
    if (!stats.is_monomorphic)
    {
        variant.array() /= stats.stddev;
    }
    return stats;
}

auto OrthStandardizeMethod::process_dominant(
    Eigen::Ref<Eigen::VectorXd> variant) -> VariantStats
{
    auto stats = OrthCenterMethod::process_dominant(variant);
    if (!stats.is_monomorphic)
    {
        variant.array() /= stats.stddev;
    }
    return stats;
}

}  // namespace gelex
