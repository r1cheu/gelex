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

#include "vi_command.h"

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
#include "gelex/data/pipe/geno.h"
#include "gelex/data/pipe/pheno.h"
#include "gelex/engine/vi.h"
#include "gelex/infra/logging/dataset_event.h"
#include "gelex/infra/logging/notify.h"
#include "gelex/model/bayes/builder.h"
#include "gelex/model/bayes/method.h"
#include "gelex/model/bayes/model.h"
#include "vi_config.h"
#include "vi_reporter.h"

auto vi_execute(argparse::ArgumentParser& cmd) -> int
{
    auto engine_config = gelex::cli::make_vi_config(cmd);
    auto [pheno_config, geno_config]
        = gelex::cli::make_fit_data_configs(cmd, cmd.get<bool>("--mmap"));

    const auto& method_config = engine_config.method;
    geno_config.mode = method_config.mode;

    int threads = cmd.get<int>("--threads");
    gelex::cli::ViReporter reporter;
    gelex::cli::DatasetReporter dataset_reporter;
    gelex::cli::PhenoReporter pheno_reporter;
    gelex::cli::GenoReporter geno_reporter;
    gelex::cli::setup_parallelization(threads);

    reporter.on_event(gelex::VIBannerEvent{});
    reporter.on_event(
        gelex::VIConfigEvent{
            .method = method_config,
            .mode = method_config.mode,
            .max_iters = static_cast<int>(engine_config.params.max_iters),
            .tol = engine_config.params.tol,
        });

    gelex::notify(dataset_reporter.as_observer(), gelex::DatasetSectionEvent{});

    gelex::PhenoPipe pheno(pheno_config, pheno_reporter.as_observer());
    pheno.load();

    gelex::GenoPipe geno(geno_config, geno_reporter.as_observer());

    std::vector<const gelex::dataframe::Index<std::string>*> all_indices{
        &geno.sample_index(), &pheno.pheno_index()};
    all_indices.append_range(pheno.covar_indices());
    auto common = gelex::cli::intersect_or_throw(
        std::move(all_indices),
        dataset_reporter.as_observer(),
        "phenotype, genotype (.fam), and covariates");

    pheno.gather(common);
    geno.load(common);

    auto model = gelex::build_bayes_model(std::move(pheno), std::move(geno));
    auto stats = gelex::compute_genetic_stats(model, method_config);
    auto method = gelex::bayes::build_bayes_method(
        method_config, stats, model.phenotype_variance());

    gelex::vi::Engine engine(std::move(engine_config));
    engine.run(model, std::move(method), reporter.as_observer());

    return 0;
}
