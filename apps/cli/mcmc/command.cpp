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

#include <utility>

#include <CLI/CLI.hpp>

#include "cli/bayes_recipe_options.h"
#include "cli/cli_helper.h"
#include "cli/data_pipe_config.h"
#include "cli/dataset_reporter.h"
#include "cli/geno_reporter.h"
#include "cli/pheno_reporter.h"
#include "config.h"
#include "gelex/bayes/builder.h"
#include "gelex/bayes/recipe.h"
#include "gelex/data/pipe/geno.h"
#include "gelex/data/pipe/pheno.h"
#include "gelex/engine/mcmc.h"
#include "gelex/infra/logging/dataset_event.h"
#include "gelex/infra/logging/fit_event.h"
#include "gelex/infra/logging/notify.h"
#include "reporter.h"

auto mcmc_execute(CLI::App& cmd) -> int
{
    auto recipe_options = cli::make_bayes_recipe_options(cmd);
    auto engine_config = cli::make_mcmc_engine_config(cmd);
    auto [pheno_config, geno_config]
        = cli::make_dataset_configs(cmd, cmd.get_option("--mmap")->count() > 0);

    int threads = cmd.get_option("--threads")->as<int>();
    cli::McmcReporter reporter;
    cli::DatasetReporter dataset_reporter;
    cli::PhenoReporter pheno_reporter;
    cli::GenoReporter geno_reporter;
    cli::setup_parallelization(threads);

    gelex::notify(reporter.as_observer(), gelex::MCMCBannerEvent{});
    gelex::notify(
        reporter.as_observer(),
        gelex::MCMCConfigEvent{
            .recipe_scheme = recipe_options.scheme,
            .requested_effects = recipe_options.modes,
            .n_iters = engine_config.mcmc_params.n_iters,
            .n_burn_in = engine_config.mcmc_params.n_burn_in,
            .seed = engine_config.seed,
        });

    gelex::notify(dataset_reporter.as_observer(), gelex::DatasetSectionEvent{});

    gelex::PhenoPipe pheno(pheno_config, pheno_reporter.as_observer());
    gelex::GenoPipe geno(geno_config, geno_reporter.as_observer());

    auto common = cli::intersect_or_throw(
        dataset_reporter.as_observer(),
        "phenotype, genotype (.fam), and covariates",
        geno.sample_indices(),
        pheno.sample_indices());

    pheno.load(common);
    geno.load(common);

    auto bayes_recipe = gelex::bayes::BayesRecipe(std::move(recipe_options));

    auto model = gelex::build_bayes_model(std::move(pheno), std::move(geno));
    auto prior = bayes_recipe.make_prior(model);

    gelex::mcmc::Engine engine(std::move(engine_config));
    engine.run(model, std::move(prior), reporter.as_observer());

    return 0;
}
