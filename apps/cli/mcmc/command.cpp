// Copyright 2026 RuLei Chen
// SPDX-License-Identifier: Apache-2.0

#include "command.h"

#include <Eigen/Core>
#include <cstddef>
#include <fmt/ranges.h>
#include <functional>
#include <optional>
#include <ranges>
#include <string>
#include <type_traits>
#include <utility>

#include "gelex/bayes/diagnostics.h"
#include "gelex/bayes/genetic/marker_covariate.h"
#include "gelex/bayes/genetic/marker_covariate_io.h"
#include "gelex/bayes/genotype/design.h"
#include "gelex/bayes/marker_effects.h"
#include "gelex/bayes/mcmc_runner.h"
#include "gelex/bayes/model.h"
#include "gelex/bayes/prior.h"
#include "gelex/bayes/sampling_plan.h"
#include "gelex/data/bed.h"
#include "gelex/data/snp_lut.h"
#include "gelex/data/snp_lut_io.h"
#include "gelex/genetic_mode.h"

#include "cli/common_data.h"
#include "cli/formatter.h"
#include "cli/mcmc/data.h"
#include "cli/mcmc/diagnostics_reporter.h"
#include "cli/mcmc/progress.h"
#include "cli/report_printer.h"
#include "cli/runtime.h"
#include "cli/summary.h"
#include "recipe.h"

namespace
{

constexpr double hpdi_prob = 0.95;

struct LoadedMcmcModel
{
    gelex::BayesModel model;
    std::string phenotype_name;
};

auto load_marker_covariate(const cli::McmcConfig& config, const gelex::Bed& bed)
    -> std::optional<gelex::bayes::MarkerCovariate>
{
    if (config.manno.empty())
    {
        return std::nullopt;
    }
    return gelex::bayes::make_marker_covariate(
        gelex::bayes::read_marker_annotation(config.manno), bed.bim());
}

auto load_mcmc_model(const cli::McmcConfig& config) -> LoadedMcmcModel
{
    auto bed = gelex::open_bed(config.bfile);
    const auto total_snps = static_cast<std::size_t>(bed.num_snps());
    cli::printer().block(cli::section("Genotype Processing:"));

    cli::McmcDataLoader loader(std::move(bed), config.random);
    auto base_data = cli::load_base_data(loader, config.base_data);
    auto design_data = std::move(loader).results();
    auto marker_covariate = load_marker_covariate(config, design_data.bed);

    cli::GenotypeProgress progress{total_snps};
    auto genetic = gelex::bayes::GeneticDesign{
        std::move(design_data.bed),
        config.mode,
        config.geno_method,
        std::move(marker_covariate),
        std::ref(progress)};
    progress.finish();
    auto model = gelex::BayesModel{
        std::move(base_data.phenotype),
        std::move(base_data.fixed_design),
        std::move(design_data.random),
        std::move(genetic)};

    return LoadedMcmcModel{
        .model = std::move(model),
        .phenotype_name = std::move(base_data.pheno_name)};
}

auto make_model_summary(const gelex::BayesModel& model) -> cli::Summary
{
    cli::Summary summary{"Model Summary"};
    summary.field("Fixed terms", "{}", model.fixed().column_names().size());
    const auto random_effect_names
        = model.random()
          | std::views::transform([](const auto& design)
                                  { return design.name(); });
    if (model.random().empty())
    {
        summary.field("Random effects", "None");
    }
    else
    {
        summary.field(
            "Random effects", "{}", fmt::join(random_effect_names, ", "));
    }
    if (const auto& marker_covariate = model.genetic().marker_covariate())
    {
        summary.field(
            "Marker annotation",
            "{}",
            fmt::join(marker_covariate->annotation_names().subspan(1), ", "));
    }
    for (const gelex::GeneticMode mode : model.genetic().each_mode())
    {
        const auto invalid_snps
            = model.genetic().cols()
              - static_cast<Eigen::Index>(
                  model.genetic().projection(mode).valid_indices().size());
        const std::string label = mode == gelex::GeneticMode::D
                                      ? "Dominance SNPs"
                                      : "Additive SNPs";
        if (invalid_snps == 0)
        {
            summary.field(label, "all valid");
        }
        else
        {
            summary.field(label, "{} excluded", invalid_snps);
        }
    }
    return summary;
}

template <typename Recipe>
auto run_mcmc(const cli::McmcConfig& config, const Recipe& recipe) -> int
{
    gelex::MCMCRunner runner{gelex::SamplingPlan{
        config.iters, config.burn_in, config.thin, config.seed}};
    auto loaded = load_mcmc_model(config);
    auto& model = loaded.model;

    gelex::ModeMap<gelex::SnpLutMatrix> snp_luts;
    for (const auto mode : model.genetic().each_mode())
    {
        snp_luts.emplace(mode, model.genetic().projection(mode).snp_luts());
    }
    gelex::write_snp_luts(config.out + ".snplut", snp_luts);

    cli::Summary{"Dataset Summary"}
        .field("Trait", "{}", loaded.phenotype_name)
        .field("Samples", "{}", model.num_individuals())
        .field("Variants", "{}", model.genetic().cols())
        .show();

    make_model_summary(model).show();

    const auto prior = gelex::make_prior(recipe, model);
    cli::printer().block(cli::section("MCMC Sampling:"));
    const auto total_iterations = static_cast<std::size_t>(config.iters);
    cli::McmcProgress progress{total_iterations, config.burn_in};
    const auto draws_path = config.out + ".draws";
    runner.run(model, prior, draws_path, std::ref(progress));
    progress.finish();

    using prior_type =
        typename std::remove_cvref_t<decltype(prior)>::genetic_prior_type;
    const auto diagnostics
        = gelex::read_diagnostics<prior_type>(draws_path, model, hpdi_prob);
    const auto entries = gelex::diagnostic_entries(diagnostics);
    gelex::write_diagnostics(config.out + ".summary", entries);
    cli::show_diagnostics(entries);
    const auto marker_effects
        = gelex::read_marker_effects<prior_type>(draws_path, model);
    gelex::write_marker_effects(
        config.out + ".snpeff", model.genetic(), marker_effects);

    cli::printer().block(
        cli::results_saved(
            config.out,
            ".draws, .draws.csc, .snplut, .summary, .snpeff, .log"));
    return 0;
}

}  // namespace

auto mcmc_execute(const cli::McmcConfig& config) -> int
{
    cli::setup_parallelization(config.threads);
    return cli::dispatch_mcmc_recipe(
        config, [&](const auto& recipe) { return run_mcmc(config, recipe); });
}
