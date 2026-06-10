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

#include "config.h"

#include <CLI/CLI.hpp>

#include <fmt/format.h>
#include <Eigen/Core>
#include <cstddef>
#include <optional>
#include <span>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

#include "cli/cli_helper.h"
#include "gelex/engine/simulation.h"
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

auto build_scheme(CLI::App& cmd, const SchemeFlags& flags, double heritability)
    -> gelex::SimulationEngine::SimulateScheme
{
    auto variances
        = cmd.get_option(flags.var_flag)->count() > 0
              ? cmd.get_option(flags.var_flag)->as<std::vector<double>>()
              : std::vector<double>{};
    auto counts = cmd.get_option(flags.n_flag)->count() > 0
                      ? cmd.get_option(flags.n_flag)->as<std::vector<int>>()
                      : std::vector<int>{};
    validate_effect_classes(variances, counts, flags.label);
    return {
        .heritability = heritability,
        .effect_sizes = create_effectsize_vec(variances, counts),
    };
}

}  // namespace

namespace cli
{

auto make_simulate_config(CLI::App& cmd) -> gelex::SimulationEngine::Config
{
    auto add_heritability
        = cmd.get_option("--h2")->count() > 0
              ? std::make_optional(cmd.get_option("--h2")->as<double>())
              : std::nullopt;
    auto dom_heritability
        = cmd.get_option("--d2")->count() > 0
              ? std::make_optional(cmd.get_option("--d2")->as<double>())
              : std::nullopt;
    auto dom_positive_prob
        = cmd.get_option("--dom-pos-prob")->count() > 0
              ? std::make_optional(
                    cmd.get_option("--dom-pos-prob")->as<double>())
              : std::nullopt;

    validate_heritabilities(
        add_heritability, dom_heritability, dom_positive_prob);

    std::optional<gelex::SimulationEngine::SimulateScheme> additive;
    if (add_heritability)
    {
        additive = build_scheme(
            cmd,
            {"--add-var", "--add-n", "Additive effect class"},
            *add_heritability);
    }

    std::optional<gelex::SimulationEngine::SimulateScheme> dominance;
    if (dom_heritability)
    {
        dominance = build_scheme(
            cmd,
            {"--dom-var", "--dom-n", "Dominance effect class"},
            *dom_heritability);
    }

    return gelex::SimulationEngine::Config{
        .seed = cmd.get_option("--seed")->as<int>(),
        .bfile_prefix = cmd.get_option("--bfile")->as<std::string>(),
        .output_prefix = cmd.get_option("--out")->as<std::string>(),
        .geno_method = parse_genotype_method(
            cmd.get_option("--geno-method")->as<std::string>()),
        .additive = std::move(additive),
        .dominance = std::move(dominance),
        .dom_positive_prob = dom_positive_prob,
    };
}

}  // namespace cli
