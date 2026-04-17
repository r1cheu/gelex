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

#include "gelex/simulate/effect_sampler.h"

#include <fmt/format.h>
#include <algorithm>
#include <cmath>
#include <optional>
#include <random>
#include <span>
#include <vector>

#include <Eigen/Core>

#include "gelex/exception.h"
#include "gelex/infra/stats/truncated_normal.h"

namespace gelex
{

EffectSampler::EffectSampler(
    std::vector<EffectSizeClass> add_classes,
    std::vector<EffectSizeClass> dom_classes,
    std::mt19937_64& rng,
    std::optional<double> dom_positive_prob)
    : add_classes_(std::move(add_classes)),
      dom_classes_(std::move(dom_classes)),
      rng_(rng),
      dom_positive_prob_(dom_positive_prob)
{
    validate_effect_classes(add_classes_, "Additive");
    if (!dom_classes_.empty())
    {
        validate_effect_classes(dom_classes_, "Dominance");
    }
}

auto EffectSampler::assign_effect_classes(
    std::span<const EffectSizeClass> classes,
    Eigen::Index n_snps) -> std::vector<int>
{
    std::vector<int> assignments(n_snps);
    const auto n_classes = static_cast<int>(classes.size());
    Eigen::Index assigned_snps = 0;

    for (int cls = 0; cls < n_classes; ++cls)
    {
        const bool is_last_cls = (cls == n_classes - 1);
        const Eigen::Index remain_snps = n_snps - assigned_snps;
        const auto expected_snps_in_cls = static_cast<Eigen::Index>(
            std::round(static_cast<double>(n_snps) * classes[cls].proportion));
        Eigen::Index class_count
            = is_last_cls ? remain_snps
                          : std::min(expected_snps_in_cls, remain_snps);

        std::fill_n(assignments.begin() + assigned_snps, class_count, cls);
        assigned_snps += class_count;
    }

    std::ranges::shuffle(assignments, rng_);
    return assignments;
}

auto EffectSampler::sample(Eigen::Index n_snps) -> CausalEffects
{
    auto add_class_assignments = assign_effect_classes(add_classes_, n_snps);
    const bool has_dominance = !dom_classes_.empty();

    auto dom_class_assignments
        = has_dominance ? assign_effect_classes(dom_classes_, n_snps)
                        : std::vector<int>{};

    CausalEffects causal_effects;
    causal_effects.resize(n_snps);
    for (Eigen::Index i = 0; i < n_snps; ++i)
    {
        const int add_class = add_class_assignments[i];
        causal_effects.additive(i)
            = sample_effect_value(add_classes_, add_class);
        causal_effects.add_class(i) = add_class;

        if (has_dominance)
        {
            const int dom_class = dom_class_assignments[i];
            causal_effects.dominance(i)
                = dom_positive_prob_.has_value()
                      ? sample_truncated_effect_value(dom_classes_, dom_class)
                      : sample_effect_value(dom_classes_, dom_class);
            causal_effects.dom_class(i) = dom_class;
        }
    }
    return causal_effects;
}

auto EffectSampler::sample_effect_value(
    std::span<const EffectSizeClass> classes,
    int cls) -> double
{
    const double variance = classes[cls].variance;
    if (variance == 0.0)
    {
        return 0.0;
    }

    auto distribution
        = std::normal_distribution<double>(0.0, std::sqrt(variance));
    return distribution(rng_);
}

auto EffectSampler::sample_truncated_effect_value(
    std::span<const EffectSizeClass> classes,
    int cls) -> double
{
    const double variance = classes[cls].variance;
    if (variance == 0.0)
    {
        return 0.0;
    }

    const double sigma = std::sqrt(variance);
    const double p = *dom_positive_prob_;
    std::bernoulli_distribution coin(p);

    if (coin(rng_))
    {
        return detail::sample_left_truncated_normal(0.0, sigma, 0.0, rng_);
    }
    return detail::sample_right_truncated_normal(0.0, sigma, 0.0, rng_);
}

void EffectSampler::validate_effect_classes(
    std::span<const EffectSizeClass> classes,
    std::string_view label)
{
    if (classes.empty())
    {
        throw GelexException(
            fmt::format("{} effect classes must not be empty", label));
    }

    double total_proportion = 0.0;
    for (const auto& cls : classes)
    {
        if (cls.proportion <= 0.0)
        {
            throw GelexException(
                fmt::format("{} effect class proportion must be > 0", label));
        }
        if (cls.variance < 0.0)
        {
            throw GelexException(
                fmt::format("{} effect class variance must be >= 0", label));
        }
        total_proportion += cls.proportion;
    }

    if (std::abs(total_proportion - 1.0) > 1e-6)
    {
        throw GelexException(
            fmt::format(
                "{} effect class proportions must sum to 1.0 (got {})",
                label,
                total_proportion));
    }
}
}  // namespace gelex
