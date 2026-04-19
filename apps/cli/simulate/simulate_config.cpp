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

#include <fmt/format.h>
#include <optional>
#include <span>
#include <string_view>
#include <utility>
#include <vector>

#include "cli/cli_helper.h"
#include "gelex/exception.h"
#include "gelex/types/sim_types.h"

namespace
{

auto validate_effect_classes(
    std::span<const double> variances,
    std::span<const int> counts,
    std::string_view label) -> void
{
    if (variances.size() != counts.size())
    {
        throw gelex::GelexException(
            fmt::format(
                "{} variances and counts must have the same number of values",
                label));
    }
}

auto create_effectsize_vec(
    std::span<const double> variances,
    std::span<const int> counts) -> std::vector<gelex::EffectSize>
{
    std::vector<gelex::EffectSize> classes(variances.size());
    for (size_t i = 0; i < variances.size(); ++i)
    {
        classes[i] = {static_cast<Eigen::Index>(counts[i]), variances[i]};
    }
    return classes;
}

}  // namespace

namespace gelex::cli
{

auto make_simulate_config(argparse::ArgumentParser& cmd)
    -> SimulationEngine::Config
{
    auto dom_heritability = cmd.is_used("--d2")
                                ? std::make_optional(cmd.get<double>("--d2"))
                                : std::nullopt;
    auto dom_positive_prob
        = cmd.is_used("--dom-pos-prob")
              ? std::make_optional(cmd.get<double>("--dom-pos-prob"))
              : std::nullopt;

    auto add_variances = cmd.get<std::vector<double>>("--add-var");
    auto add_counts = cmd.get<std::vector<int>>("--add-n");
    auto dom_variances = cmd.get<std::vector<double>>("--dom-var");
    auto dom_counts = cmd.is_used("--dom-n")
                          ? cmd.get<std::vector<int>>("--dom-n")
                          : std::vector<int>{};

    auto add_heritability = cmd.get<double>("--h2");

    if (add_heritability <= 0.0 || add_heritability >= 1.0)
    {
        throw gelex::GelexException("Heritability must be in (0, 1)");
    }
    if (dom_heritability
        && (*dom_heritability < 0.0 || *dom_heritability >= 1.0))
    {
        throw gelex::GelexException(
            "Dominance variance (d2) must be in [0, 1)");
    }
    if (dom_heritability && add_heritability + *dom_heritability >= 1.0)
    {
        throw gelex::GelexException("h2 + d2 must be less than 1");
    }
    if (dom_positive_prob
        && (*dom_positive_prob <= 0.0 || *dom_positive_prob >= 1.0))
    {
        throw gelex::GelexException(
            "Dominance positive probability must be in (0, 1)");
    }
    if (dom_positive_prob && !dom_heritability)
    {
        throw gelex::GelexException(
            "--dom-pos-prob requires --d2 to be specified");
    }

    validate_effect_classes(add_variances, add_counts, "Additive effect class");
    if (dom_heritability)
    {
        validate_effect_classes(
            dom_variances, dom_counts, "Dominance effect class");
    }

    std::optional<SimulationEngine::SimulateScheme> dominance;
    if (dom_heritability)
    {
        dominance = SimulationEngine::SimulateScheme{
            .heritability = *dom_heritability,
            .effect_sizes = create_effectsize_vec(dom_variances, dom_counts),
        };
    }

    return SimulationEngine::Config{
        .seed = cmd.get<int>("--seed"),
        .bfile_prefix = cmd.get("--bfile"),
        .output_prefix = cmd.get("--out"),
        .geno_method
        = parse_genotype_process_method(cmd.get<std::string>("--geno-method")),
        .additive
        = SimulationEngine::SimulateScheme{
            .heritability = add_heritability,
            .effect_sizes = create_effectsize_vec(add_variances, add_counts),
        },
        .dominance = std::move(dominance),
        .dom_positive_prob = dom_positive_prob,
    };
}

}  // namespace gelex::cli
