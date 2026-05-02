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
#include "gelex/simulate/sim_types.h"

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

auto validate_fraction_open(std::optional<double> value, std::string_view label)
    -> void
{
    if (value && (*value <= 0.0 || *value >= 1.0))
    {
        throw gelex::GelexException(fmt::format("{} must be in (0, 1)", label));
    }
}

auto validate_heritabilities(
    std::optional<double> add_h2,
    std::optional<double> dom_h2,
    std::optional<double> dom_pos_prob) -> void
{
    if (!add_h2 && !dom_h2)
    {
        throw gelex::GelexException(
            "Must specify at least one of --h2 or --d2");
    }
    validate_fraction_open(add_h2, "Heritability (h2)");
    validate_fraction_open(dom_h2, "Dominance variance (d2)");
    if (add_h2 && dom_h2 && *add_h2 + *dom_h2 >= 1.0)
    {
        throw gelex::GelexException("h2 + d2 must be less than 1");
    }
    validate_fraction_open(dom_pos_prob, "Dominance positive probability");
    if (dom_pos_prob && !dom_h2)
    {
        throw gelex::GelexException(
            "--dom-pos-prob requires --d2 to be specified");
    }
}

struct SchemeFlags
{
    std::string var_flag;
    std::string n_flag;
    std::string label;
};

auto build_scheme(
    argparse::ArgumentParser& cmd,
    const SchemeFlags& flags,
    double heritability) -> gelex::SimulationEngine::SimulateScheme
{
    auto variances = cmd.is_used(flags.var_flag)
                         ? cmd.get<std::vector<double>>(flags.var_flag)
                         : std::vector<double>{};
    auto counts = cmd.is_used(flags.n_flag)
                      ? cmd.get<std::vector<int>>(flags.n_flag)
                      : std::vector<int>{};
    validate_effect_classes(variances, counts, flags.label);
    return {
        .heritability = heritability,
        .effect_sizes = create_effectsize_vec(variances, counts),
    };
}

}  // namespace

namespace gelex::cli
{

auto make_simulate_config(argparse::ArgumentParser& cmd)
    -> SimulationEngine::Config
{
    auto add_heritability = cmd.is_used("--h2")
                                ? std::make_optional(cmd.get<double>("--h2"))
                                : std::nullopt;
    auto dom_heritability = cmd.is_used("--d2")
                                ? std::make_optional(cmd.get<double>("--d2"))
                                : std::nullopt;
    auto dom_positive_prob
        = cmd.is_used("--dom-pos-prob")
              ? std::make_optional(cmd.get<double>("--dom-pos-prob"))
              : std::nullopt;

    validate_heritabilities(
        add_heritability, dom_heritability, dom_positive_prob);

    std::optional<SimulationEngine::SimulateScheme> additive;
    if (add_heritability)
    {
        additive = build_scheme(
            cmd,
            {"--add-var", "--add-n", "Additive effect class"},
            *add_heritability);
    }

    std::optional<SimulationEngine::SimulateScheme> dominance;
    if (dom_heritability)
    {
        dominance = build_scheme(
            cmd,
            {"--dom-var", "--dom-n", "Dominance effect class"},
            *dom_heritability);
    }

    return SimulationEngine::Config{
        .seed = cmd.get<int>("--seed"),
        .bfile_prefix = cmd.get("--bfile"),
        .output_prefix = cmd.get("--out"),
        .geno_method
        = parse_genotype_process_method(cmd.get<std::string>("--geno-method")),
        .additive = std::move(additive),
        .dominance = std::move(dominance),
        .dom_positive_prob = dom_positive_prob,
    };
}

}  // namespace gelex::cli
