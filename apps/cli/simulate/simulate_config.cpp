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

}  // namespace

namespace gelex::cli
{

auto make_simulate_config(argparse::ArgumentParser& cmd) -> SimulateConfig
{
    auto dom_heritability = cmd.is_used("--d2")
                                ? std::make_optional(cmd.get<double>("--d2"))
                                : std::nullopt;

    SimulateConfig config{
        .bed_path = gelex::format_bed_path(cmd.get("--bfile")),
        .output_path = cmd.get("--out"),

        .intercept = cmd.get<double>("--intercept"),

        .add_heritability = cmd.get<double>("--h2"),
        .additive_variances = cmd.get<std::vector<double>>("--add-var"),
        .additive_proportions = cmd.get<std::vector<double>>("--add-prop"),

        .dom_heritability = dom_heritability,
        .dominance_variances = cmd.get<std::vector<double>>("--dom-var"),
        .dominance_proportions = cmd.get<std::vector<double>>("--dom-prop"),

        .seed = cmd.get<int>("--seed")};

    if (config.add_heritability <= 0.0 || config.add_heritability >= 1.0)
    {
        throw gelex::ArgumentValidationException(
            "Heritability must be in (0, 1)");
    }
    if (config.dom_heritability
        && (*config.dom_heritability < 0.0 || *config.dom_heritability >= 1.0))
    {
        throw gelex::ArgumentValidationException(
            "Dominance variance (d2) must be in [0, 1)");
    }
    if (config.dom_heritability
        && config.add_heritability + *config.dom_heritability >= 1.0)
    {
        throw gelex::ArgumentValidationException("h2 + d2 must be less than 1");
    }

    validate_effect_classes(
        config.additive_variances,
        config.additive_proportions,
        "Additive effect class");
    if (config.dom_heritability)
    {
        validate_effect_classes(
            config.dominance_variances,
            config.dominance_proportions,
            "Dominance effect class");
    }

    return config;
}

}  // namespace gelex::cli
