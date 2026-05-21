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

#include "mcmc_command.h"

#include <string>
#include <utility>

#include <argparse.h>

#include "cli/bayes_recipe_config.h"
#include "cli/cli_helper.h"
#include "cli/data_pipe_config.h"
#include "cli/dataset_reporter.h"
#include "cli/geno_reporter.h"
#include "cli/pheno_reporter.h"
#include "gelex/data/pipe/geno.h"
#include "gelex/data/pipe/pheno.h"
#include "gelex/engine/mcmc.h"
#include "gelex/infra/logging/dataset_event.h"
#include "gelex/infra/logging/notify.h"
#include "gelex/model/bayes/builder.h"
#include "gelex/model/bayes/legacy_method.h"
#include "gelex/model/bayes/model.h"
#include "gelex/model/bayes/recipe.h"
#include "mcmc_config.h"
#include "mcmc_overrides.h"
#include "mcmc_reporter.h"

auto mcmc_execute(argparse::ArgumentParser& cmd) -> int
{
    auto mcmc_config = gelex::cli::make_mcmc_config(cmd);
    auto [pheno_config, geno_config]
        = gelex::cli::make_dataset_configs(cmd, cmd.get<bool>("--mmap"));

    int threads = cmd.get<int>("--threads");
    gelex::cli::McmcReporter reporter;
    gelex::cli::DatasetReporter dataset_reporter;
    gelex::cli::PhenoReporter pheno_reporter;
    gelex::cli::GenoReporter geno_reporter;
    gelex::cli::setup_parallelization(threads);

    gelex::notify(reporter.as_observer(), gelex::MCMCBannerEvent{});
    gelex::notify(
        reporter.as_observer(),
        gelex::MCMCConfigEvent{
            .method = mcmc_config.engine.method,
            .requested_effects = mcmc_config.engine.requested_effects,
            .n_iters = mcmc_config.engine.mcmc_params.n_iters,
            .n_burn_in = mcmc_config.engine.mcmc_params.n_burn_in,
            .seed = mcmc_config.engine.seed,
        });

    gelex::notify(dataset_reporter.as_observer(), gelex::DatasetSectionEvent{});

    gelex::PhenoPipe pheno(pheno_config, pheno_reporter.as_observer());
    gelex::GenoPipe geno(geno_config, geno_reporter.as_observer());

    auto common = gelex::cli::intersect_or_throw(
        dataset_reporter.as_observer(),
        "phenotype, genotype (.fam), and covariates",
        geno.sample_indices(),
        pheno.sample_indices());

    pheno.load(common);
    geno.load(common);

    auto bayes_config = gelex::cli::make_bayes_recipe_config(cmd);
    auto preset = gelex::bayes::to_bayes_recipe_preset(
        cmd.get<std::string>("--method"));
    auto bayes_recipe = gelex::bayes::BayesRecipe(preset, bayes_config);

    auto model = gelex::build_bayes_model(std::move(pheno), std::move(geno));
    auto stats = gelex::compute_genetic_stats(model);
    auto prior = bayes_recipe.make_prior(model);

    auto method = gelex::bayes::build_bayes_method(
        mcmc_config.engine.method, stats, model.phenotype_variance());
    gelex::cli::apply_overrides(method, mcmc_config.overrides);

    gelex::mcmc::Engine engine(std::move(mcmc_config.engine));
    engine.run(model, std::move(method), reporter.as_observer());

    return 0;
}
