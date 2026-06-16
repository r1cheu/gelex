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

#include "command.h"

#include <algorithm>
#include <cstddef>
#include <optional>
#include <random>
#include <span>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

#include <fmt/format.h>
#include <CLI/CLI.hpp>
#include <Eigen/Core>

#include "cli/cli_helper.h"
#include "gelex/data/reader.h"
#include "gelex/data/sample_id.h"
#include "gelex/exception.h"
#include "gelex/infra/logging/notify.h"
#include "gelex/infra/logging/simulate_event.h"
#include "gelex/infra/stats/detail/var.h"
#include "gelex/io/detail/text_writer.h"
#include "gelex/simulate/effect_sampler.h"
#include "gelex/simulate/genetic_value_calculator.h"
#include "gelex/simulate/genetic_value_scaler.h"
#include "gelex/simulate/sim_types.h"
#include "gelex/types/genetic_effect_type.h"
#include "reporter.h"

namespace
{

struct SimulateScheme
{
    double heritability;
    std::vector<gelex::EffectSize> effect_sizes;
};

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
    for (std::size_t i = 0; i < variances.size(); ++i)
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
    -> SimulateScheme
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

auto simulate_execute(CLI::App& sim) -> int
{
    auto add_heritability
        = sim.get_option("--h2")->count() > 0
              ? std::make_optional(sim.get_option("--h2")->as<double>())
              : std::nullopt;
    auto dom_heritability
        = sim.get_option("--d2")->count() > 0
              ? std::make_optional(sim.get_option("--d2")->as<double>())
              : std::nullopt;
    auto dom_positive_prob
        = sim.get_option("--dom-pos-prob")->count() > 0
              ? std::make_optional(
                    sim.get_option("--dom-pos-prob")->as<double>())
              : std::nullopt;

    validate_heritabilities(
        add_heritability, dom_heritability, dom_positive_prob);

    std::optional<SimulateScheme> additive_scheme;
    if (add_heritability)
    {
        additive_scheme = build_scheme(
            sim,
            {"--add-var", "--add-n", "Additive effect class"},
            *add_heritability);
    }

    std::optional<SimulateScheme> dominance_scheme;
    if (dom_heritability)
    {
        dominance_scheme = build_scheme(
            sim,
            {"--dom-var", "--dom-n", "Dominance effect class"},
            *dom_heritability);
    }

    cli::SimulatorReporter reporter;
    auto observer = reporter.as_observer();

    gelex::notify(observer, gelex::SimulateBannerEvent{});
    gelex::notify(
        observer,
        gelex::SimulateConfigLoadedEvent{
            .add_heritability = add_heritability,
            .dom_heritability = dom_heritability,
            .seed = sim.get_option("--seed")->as<int>(),
        });

    std::mt19937_64 rng(sim.get_option("--seed")->as<int>());

    const auto bfile_prefix = sim.get_option("--bfile")->as<std::string>();
    auto bim = gelex::read_bim(bfile_prefix + ".bim");
    auto fam = gelex::read_fam(bfile_prefix + ".fam");

    gelex::GeneticValueCalculator calculator(bfile_prefix, bim, fam);

    auto all_ids = bim.index().keys();
    std::vector<std::string_view> shuffled_ids(all_ids.begin(), all_ids.end());
    std::ranges::shuffle(shuffled_ids, rng);

    const auto geno_method = cli::parse_genotype_method(
        sim.get_option("--geno-method")->as<std::string>());

    std::optional<gelex::GeneticValues> additive;
    if (additive_scheme)
    {
        additive.emplace();
        additive->causal_snps = gelex::NormalSampler(
            shuffled_ids, additive_scheme->effect_sizes)(rng);
        additive->gebv = calculator.calculate<gelex::GeneticMode::A>(
            *additive, geno_method, observer);
    }

    std::optional<gelex::GeneticValues> dominance;
    if (dominance_scheme)
    {
        dominance.emplace();
        dominance->causal_snps = gelex::NormalSampler(
            shuffled_ids, dominance_scheme->effect_sizes)(rng);
        dominance->gebv = calculator.calculate<gelex::GeneticMode::D>(
            *dominance, geno_method, observer);
    }

    const Eigen::Index n_samples
        = additive ? additive->gebv.size() : dominance->gebv.size();
    std::normal_distribution<double> normal01(0.0, 1.0);
    Eigen::VectorXd residual = Eigen::VectorXd::NullaryExpr(
        n_samples, [&] { return normal01(rng); });

    auto h2_of = [](const auto& scheme) -> std::optional<double>
    {
        return scheme ? std::optional<double>(scheme->heritability)
                      : std::nullopt;
    };
    gelex::GeneticValueScaler scaler(
        h2_of(additive_scheme), h2_of(dominance_scheme));
    scaler.scale(
        additive ? &*additive : nullptr,
        dominance ? &*dominance : nullptr,
        residual);

    Eigen::VectorXd phenotypes = additive ? additive->gebv : dominance->gebv;
    if (additive && dominance)
    {
        phenotypes += dominance->gebv;
    }
    phenotypes += residual;

    const double var_phen = gelex::stats::detail::vecvar(
        phenotypes, gelex::stats::detail::VarNormType::Population);
    gelex::SimulateVarianceSummaryEvent summary;
    if (additive && var_phen > 0.0)
    {
        summary.realized_h2
            = gelex::stats::detail::vecvar(
                  additive->gebv, gelex::stats::detail::VarNormType::Population)
              / var_phen;
    }
    if (dominance && var_phen > 0.0)
    {
        summary.realized_d2 = gelex::stats::detail::vecvar(
                                  dominance->gebv,
                                  gelex::stats::detail::VarNormType::Population)
                              / var_phen;
    }
    gelex::notify(observer, summary);

    const auto out_prefix = sim.get_option("--out")->as<std::string>();
    gelex::io::detail::TextWriter writer(out_prefix + ".phen");
    writer.write_header({"FID", "IID", "Phenotype"});
    for (Eigen::Index i = 0; i < n_samples; ++i)
    {
        auto [fid, iid] = gelex::split_sample_id(fam.index().keys()[i]);
        writer.write(fmt::format("{}\t{}\t{}", fid, iid, phenotypes(i)));
    }

    gelex::io::detail::TextWriter effect_writer(out_prefix + ".causal");
    auto write_single =
        [&](std::string_view column, const std::vector<gelex::CausalSnp>& snps)
    {
        effect_writer.write_header({"id", std::string(column)});
        for (const auto& snp : snps)
        {
            effect_writer.write(fmt::format("{}\t{}", snp.id, snp.effect));
        }
    };
    if (additive && dominance)
    {
        effect_writer.write_header({"id", "additive", "dominance"});
        const auto& add_snps = additive->causal_snps;
        const auto& dom_snps = dominance->causal_snps;
        for (std::size_t i = 0; i < add_snps.size(); ++i)
        {
            effect_writer.write(
                fmt::format(
                    "{}\t{}\t{}",
                    add_snps[i].id,
                    add_snps[i].effect,
                    dom_snps[i].effect));
        }
    }
    else if (additive)
    {
        write_single("additive", additive->causal_snps);
    }
    else
    {
        write_single("dominance", dominance->causal_snps);
    }

    return 0;
}
