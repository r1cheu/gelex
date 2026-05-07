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

#include <filesystem>
#include <string>
#include <utility>
#include <vector>

#include <argparse.h>

#include "cli/cli_helper.h"
#include "cli/data_pipe_config.h"
#include "cli/dataset_loader.h"
#include "cli/dataset_reporter.h"
#include "cli/geno_reporter.h"
#include "cli/pheno_reporter.h"
#include "gelex/data/dataframe/index.h"
#include "gelex/data/genotype/bed_path.h"
#include "gelex/data/pipe/geno.h"
#include "gelex/data/pipe/pheno.h"
#include "gelex/data/reader.h"
#include "gelex/engine/mcmc.h"
#include "gelex/infra/logging/dataset_event.h"
#include "gelex/infra/logging/notify.h"
#include "gelex/model/bayes/builder.h"
#include "gelex/model/bayes/method.h"
#include "gelex/model/bayes/model.h"
#include "mcmc_config.h"
#include "mcmc_overrides.h"
#include "mcmc_reporter.h"

auto mcmc_execute(argparse::ArgumentParser& cmd) -> int
{
    auto mcmc_config = gelex::cli::make_mcmc_config(cmd);
    auto [pheno_config, geno_config]
        = gelex::cli::make_fit_data_configs(cmd, cmd.get<bool>("--mmap"));

    const auto& method_config = mcmc_config.engine.method;
    geno_config.mode = method_config.mode;

    int threads = cmd.get<int>("--threads");
    gelex::cli::McmcReporter reporter;
    gelex::cli::DatasetReporter dataset_reporter;
    gelex::cli::PhenoReporter pheno_reporter;
    gelex::cli::GenoReporter geno_reporter;
    gelex::cli::setup_parallelization(threads);

    reporter.on_event(gelex::MCMCBannerEvent{});
    reporter.on_event(
        gelex::MCMCConfigEvent{
            .method = method_config,
            .mode = method_config.mode,
            .n_iters = static_cast<int>(mcmc_config.engine.mcmc_params.n_iters),
            .n_burn_in
            = static_cast<int>(mcmc_config.engine.mcmc_params.n_burn_in),
            .seed = mcmc_config.engine.seed,
        });

    gelex::notify(dataset_reporter.as_observer(), gelex::DatasetSectionEvent{});

    auto bed_path
        = gelex::genotype::format_bed_path(cmd.get<std::string>("--bfile"));
    auto fam_index
        = gelex::read_fam(
              std::filesystem::path(bed_path).replace_extension(".fam"))
              .index();

    gelex::PhenoPipe pheno(pheno_config, pheno_reporter.as_observer());
    pheno.load();

    std::vector<const gelex::dataframe::Index<std::string>*> all_indices{
        &fam_index, &pheno.pheno_index()};
    all_indices.append_range(pheno.covar_indices());
    auto common = gelex::cli::intersect_or_throw(
        std::move(all_indices),
        dataset_reporter.as_observer(),
        "phenotype, genotype (.fam), and covariates");

    pheno.gather(common);

    gelex::GenoPipe geno(geno_config, geno_reporter.as_observer());
    geno.load(common);

    auto model = gelex::build_bayes_model(std::move(pheno), std::move(geno));
    auto stats = gelex::compute_genetic_stats(model, method_config);
    auto method = gelex::bayes::build_bayes_method(
        method_config, stats, model.phenotype_variance());
    gelex::cli::apply_overrides(method, mcmc_config.overrides);

    gelex::mcmc::Engine engine(std::move(mcmc_config.engine));
    engine.run(model, std::move(method), reporter.as_observer());

    return 0;
}
