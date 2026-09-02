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

#include <string>
#include <utility>

#include "gelex/bayes/draws.h"
#include "gelex/bayes/genotype/design.h"
#include "gelex/bayes/genotype/progress.h"
#include "gelex/bayes/mcmc_progress.h"
#include "gelex/bayes/mcmc_runner.h"
#include "gelex/bayes/model.h"
#include "gelex/bayes/prior.h"
#include "gelex/bayes/result.h"
#include "gelex/bayes/result_io.h"
#include "gelex/data/bed.h"
#include "gelex/data/snp_lut.h"
#include "gelex/data/snp_lut_io.h"
#include "gelex/genetic_mode.h"

#include "cli/common_data.h"
#include "cli/formatter.h"
#include "cli/mcmc/data.h"
#include "cli/report_printer.h"
#include "cli/runtime.h"
#include "recipe.h"
#include "reporter.h"

namespace
{

struct LoadedMcmcModel
{
    gelex::BayesModel model;
    std::string phenotype_name;
};

auto load_mcmc_model(const cli::McmcConfig& config, cli::GenoReporter& reporter)
    -> LoadedMcmcModel
{
    auto bed = gelex::open_bed(config.bfile);
    cli::printer().block(cli::section("Genotype Processing:"));
    reporter.show_total(bed.num_snps());

    cli::McmcDataLoader loader(std::move(bed), config.random);
    auto base_data = cli::load_base_data(loader, config.base_data);
    auto design_data = std::move(loader).results();

    const gelex::GenoObserver observer
        = [&reporter](const gelex::GenotypeProgressEvent& event)
    { reporter.on_event(event); };
    auto genetic = gelex::bayes::GeneticDesign{
        std::move(design_data.bed), config.mode, config.geno_method, observer};
    auto model = gelex::BayesModel{
        std::move(base_data.phenotype),
        std::move(base_data.fixed_design),
        std::move(design_data.random),
        std::move(genetic)};
    reporter.show_loaded(model.genetic());

    return LoadedMcmcModel{
        .model = std::move(model),
        .phenotype_name = std::move(base_data.pheno_name)};
}

template <typename Recipe>
auto run_mcmc(
    const cli::McmcConfig& config,
    const Recipe& recipe,
    cli::GenoReporter& geno_reporter,
    cli::McmcReporter& mcmc_reporter) -> int
{
    gelex::MCMCRunner runner{config.iters, config.burn_in, config.thin};
    auto loaded = load_mcmc_model(config, geno_reporter);
    auto& model = loaded.model;
    mcmc_reporter.show_dataset_summary(model, loaded.phenotype_name);

    const auto prior = gelex::make_prior(recipe, model);
    const gelex::MCMCObserver observer
        = [&mcmc_reporter](const gelex::MCMCProgressEvent& event)
    { mcmc_reporter.on_event(event); };
    const auto result = [&]()
    {
        gelex::BayesDraws draws(
            prior, model, config.out + ".draws", runner.draw_count());
        runner.run(model, prior, draws, config.seed, observer);
        return gelex::make_result(model, draws);
    }();

    gelex::write_params(result, config.out);
    gelex::write_summary(result, config.out);

    gelex::ModeMap<gelex::SnpLutMatrix> snp_luts;
    for (const auto mode : model.genetic().each_mode())
    {
        snp_luts.emplace(mode, model.genetic().projection(mode).snp_luts());
    }
    gelex::write_snp_luts(config.out + ".snplut", snp_luts);

    // TODO: Restore .snpeff after typed SNP-effect serialization is defined.
    cli::printer().block(
        cli::results_saved(
            config.out, ".draws, .params, .summary, .snplut, .log"));
    return 0;
}

}  // namespace

auto mcmc_execute(const cli::McmcConfig& config) -> int
{
    cli::setup_parallelization(config.threads);
    cli::GenoReporter geno_reporter;
    cli::McmcReporter mcmc_reporter;
    return cli::dispatch_mcmc_recipe(
        config,
        [&](const auto& recipe)
        { return run_mcmc(config, recipe, geno_reporter, mcmc_reporter); });
}
