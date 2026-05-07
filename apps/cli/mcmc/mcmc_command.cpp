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
#include "cli/data_pipe_reporter.h"
#include "cli/fit_reporter.h"
#include "gelex/data/dataframe/index.h"
#include "gelex/data/genotype/bed_path.h"
#include "gelex/data/pipe/geno.h"
#include "gelex/data/pipe/pheno.h"
#include "gelex/data/reader.h"
#include "gelex/engine/mcmc.h"
#include "gelex/exception.h"
#include "gelex/infra/logging/data_pipe_event.h"
#include "gelex/infra/logging/notify.h"
#include "gelex/model/bayes/builder.h"
#include "gelex/model/bayes/method.h"
#include "gelex/model/bayes/model.h"
#include "mcmc_config.h"
#include "mcmc_overrides.h"

auto mcmc_execute(argparse::ArgumentParser& cmd) -> int
{
    auto fit_config = gelex::cli::make_mcmc_fit_config(cmd);
    auto [pheno_config, geno_config]
        = gelex::cli::make_fit_data_configs(cmd, cmd.get<bool>("--mmap"));

    const auto& method_config = fit_config.engine.method;
    geno_config.mode = method_config.mode;

    int threads = cmd.get<int>("--threads");
    gelex::cli::FitReporter reporter;
    gelex::cli::DataPipeReporter data_reporter;
    gelex::cli::setup_parallelization(threads);

    reporter.on_event(gelex::FitMCMCBannerEvent{});
    reporter.on_event(
        gelex::FitMCMCConfigEvent{
            .method = method_config,
            .mode = method_config.mode,
            .n_iters = static_cast<int>(fit_config.engine.mcmc_params.n_iters),
            .n_burn_in
            = static_cast<int>(fit_config.engine.mcmc_params.n_burn_in),
            .seed = fit_config.engine.seed,
        });

    auto bed_path
        = gelex::genotype::format_bed_path(cmd.get<std::string>("--bfile"));
    auto fam_index
        = gelex::read_fam(
              std::filesystem::path(bed_path).replace_extension(".fam"))
              .index();

    gelex::PhenoPipe pheno(pheno_config, data_reporter.as_observer());
    pheno.load();

    std::vector<const gelex::dataframe::Index<std::string>*> all_indices{
        &fam_index, &pheno.pheno_index()};
    all_indices.append_range(pheno.covar_indices());
    auto common = gelex::dataframe::intersect<std::string>(all_indices);

    gelex::notify(
        data_reporter.as_observer(),
        gelex::IntersectionEvent{.common_samples = common.size()});

    if (common.size() == 0)
    {
        throw gelex::GelexException(
            "No common samples across phenotype, genotype (.fam), and "
            "covariates. Check that sample IDs match across input files.");
    }

    pheno.gather(common);

    gelex::GenoPipe geno(geno_config, data_reporter.as_observer());
    geno.load(common);

    auto model = gelex::build_bayes_model(std::move(pheno), std::move(geno));
    auto stats = gelex::compute_genetic_stats(model, method_config);
    auto method = gelex::bayes::build_bayes_method(
        method_config, stats, model.phenotype_variance());
    gelex::cli::apply_overrides(method, fit_config.overrides);

    gelex::mcmc::FitEngine engine(std::move(fit_config.engine));
    engine.run(model, std::move(method), reporter.as_observer());

    return 0;
}
