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

#include "gelex/data/effect_sampler.h"

#include <algorithm>
#include <cmath>
#include <format>
#include <random>
#include <span>
#include <vector>

#include <Eigen/Core>

#include "gelex/exception.h"

namespace gelex
{

EffectSampler::EffectSampler(Config config)
    : config_(std::move(config)), rng_(config_.seed)
{
    validate_effect_classes(config_.add_classes, "Additive");
    if (config_.has_dominance)
    {
        validate_effect_classes(config_.dom_classes, "Dominance");
    }
}

auto EffectSampler::assign_effect_classes(
    std::span<const EffectSizeClass> classes,
    Eigen::Index count) -> std::vector<int>
{
    std::vector<int> assignments(count);
    const auto n_classes = static_cast<int>(classes.size());
    Eigen::Index offset = 0;

    for (int cls = 0; cls < n_classes; ++cls)
    {
        const bool is_last = (cls == n_classes - 1);
        Eigen::Index class_count
            = is_last ? count - offset
                      : std::min(
                            static_cast<Eigen::Index>(std::round(
                                static_cast<double>(count)
                                * classes[cls].proportion)),
                            count - offset);

        std::fill_n(assignments.begin() + offset, class_count, cls);
        offset += class_count;
    }

    std::ranges::shuffle(assignments, rng_);
    return assignments;
}

auto EffectSampler::sample(Eigen::Index n_snps)
    -> std::unordered_map<Eigen::Index, CausalEffect>
{
    auto add_assignments = assign_effect_classes(config_.add_classes, n_snps);

    auto dom_assignments
        = config_.has_dominance
              ? assign_effect_classes(config_.dom_classes, n_snps)
              : std::vector<int>{};

    auto sample_effect
        = [&](const std::vector<EffectSizeClass>& classes, int cls) -> double
    {
        double variance = classes[cls].variance;
        if (variance == 0.0)
        {
            return 0.0;
        }
        return std::normal_distribution<double>(0.0, std::sqrt(variance))(rng_);
    };

    std::unordered_map<Eigen::Index, CausalEffect> causal_effects;
    causal_effects.reserve(n_snps);

    for (Eigen::Index i = 0; i < n_snps; ++i)
    {
        CausalEffect effect{
            .additive = sample_effect(config_.add_classes, add_assignments[i]),
            .dominance = 0.0,
            .add_class = add_assignments[i],
            .dom_class = 0,
        };

        if (config_.has_dominance)
        {
            effect.dominance
                = sample_effect(config_.dom_classes, dom_assignments[i]);
            effect.dom_class = dom_assignments[i];
        }

        causal_effects.emplace(i, effect);
    }

    return causal_effects;
}

void EffectSampler::validate_effect_classes(
    std::span<const EffectSizeClass> classes,
    std::string_view label)
{
    if (classes.empty())
    {
        throw ArgumentValidationException(
            std::format("{} effect classes must not be empty", label));
    }

    double total_proportion = 0.0;
    for (const auto& cls : classes)
    {
        if (cls.proportion <= 0.0)
        {
            throw ArgumentValidationException(
                std::format("{} effect class proportion must be > 0", label));
        }
        if (cls.variance < 0.0)
        {
            throw ArgumentValidationException(
                std::format("{} effect class variance must be >= 0", label));
        }
        total_proportion += cls.proportion;
    }

    if (std::abs(total_proportion - 1.0) > 1e-6)
    {
        throw ArgumentValidationException(
            std::format(
                "{} effect class proportions must sum to 1.0 (got {})",
                label,
                total_proportion));
    }
}
}  // namespace gelex
