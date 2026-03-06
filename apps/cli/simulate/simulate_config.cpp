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

#include "simulate_config.h"

#include <argparse.h>

#include <format>
#include <span>
#include <string_view>
#include <vector>

#include "gelex/algo/sim/effect_sampler.h"
#include "gelex/data/genotype/bed_path.h"
#include "gelex/exception.h"

namespace
{

auto validate_effect_classes(
    std::span<const double> variances,
    std::span<const double> proportions,
    std::string_view label) -> void
{
    if (variances.size() != proportions.size())
    {
        throw gelex::ArgumentValidationException(
            std::format(
                "{} variances and proportions must have the same number of "
                "values",
                label));
    }
}

auto create_effectsize_vec(
    std::span<const double> variances,
    std::span<const double> proportions) -> std::vector<gelex::EffectSizeClass>
{
    std::vector<gelex::EffectSizeClass> classes(variances.size());
    for (size_t i = 0; i < variances.size(); ++i)
    {
        classes[i] = {proportions[i], variances[i]};
    }
    return classes;
}

}  // namespace

namespace gelex::cli
{

auto make_simulate_config(argparse::ArgumentParser& cmd)
    -> PhenotypeSimulationEngine::Config
{
    auto dom_heritability = cmd.is_used("--d2")
                                ? std::make_optional(cmd.get<double>("--d2"))
                                : std::nullopt;

    auto add_variances = cmd.get<std::vector<double>>("--add-var");
    auto add_proportions = cmd.get<std::vector<double>>("--add-prop");
    auto dom_variances = cmd.get<std::vector<double>>("--dom-var");
    auto dom_proportions = cmd.get<std::vector<double>>("--dom-prop");

    auto add_heritability = cmd.get<double>("--h2");

    if (add_heritability <= 0.0 || add_heritability >= 1.0)
    {
        throw gelex::ArgumentValidationException(
            "Heritability must be in (0, 1)");
    }
    if (dom_heritability
        && (*dom_heritability < 0.0 || *dom_heritability >= 1.0))
    {
        throw gelex::ArgumentValidationException(
            "Dominance variance (d2) must be in [0, 1)");
    }
    if (dom_heritability && add_heritability + *dom_heritability >= 1.0)
    {
        throw gelex::ArgumentValidationException("h2 + d2 must be less than 1");
    }

    validate_effect_classes(
        add_variances, add_proportions, "Additive effect class");
    if (dom_heritability)
    {
        validate_effect_classes(
            dom_variances, dom_proportions, "Dominance effect class");
    }

    return PhenotypeSimulationEngine::Config{
        .bed_path = gelex::format_bed_path(cmd.get("--bfile")),
        .output_path = cmd.get("--out"),

        .intercept = cmd.get<double>("--intercept"),

        .add_heritability = add_heritability,
        .add_effect_classes
        = create_effectsize_vec(add_variances, add_proportions),

        .dom_heritability = dom_heritability,
        .dom_effect_classes
        = dom_heritability
              ? create_effectsize_vec(dom_variances, dom_proportions)
              : std::vector<gelex::EffectSizeClass>{},

        .seed = cmd.get<int>("--seed"),
    };
}

}  // namespace gelex::cli
