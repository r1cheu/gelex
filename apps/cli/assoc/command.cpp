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

#include <array>
#include <filesystem>
#include <ranges>
#include <string>
#include <utility>
#include <vector>

#include <CLI/CLI.hpp>

#include "cli/cli_helper.h"
#include "cli/data_pipe_config.h"
#include "cli/dataset_reporter.h"
#include "cli/grm_pipe_reporter.h"
#include "cli/pheno_reporter.h"
#include "cli/reml_reporter.h"
#include "config.h"
#include "gelex/data/dataframe/index.h"
#include "gelex/data/pipe/grm.h"
#include "gelex/data/pipe/pheno.h"
#include "gelex/data/reader.h"
#include "gelex/engine/assoc.h"
#include "gelex/infra/logging/assoc_event.h"
#include "gelex/infra/logging/dataset_event.h"
#include "gelex/infra/logging/notify.h"
#include "reporter.h"

auto assoc_execute(CLI::App& cmd) -> int
{
    int threads = cmd.get_option("--threads")->as<int>();
    cli::setup_parallelization(threads);

    auto pheno_config = cli::make_pheno_config(cmd);
    pheno_config.transform_type = cli::parse_transform_type(
        cmd.get_option("--transform")->as<std::string>());
    pheno_config.int_offset = cmd.get_option("--int-offset")->as<double>();

    cli::AssocReporter reporter;
    cli::DatasetReporter dataset_reporter;
    cli::PhenoReporter pheno_reporter;
    cli::GrmPipeReporter grm_reporter;

    auto config = cli::make_assoc_config(cmd);

    reporter.on_event(gelex::AssocBannerEvent{});
    reporter.on_event(
        gelex::AssocConfigLoadedEvent{
            .mode = config.mode,
            .test_type = config.test_type,
            .loco = config.loco,
            .geno_method = config.method,
            .max_iter = config.max_iter,
            .tol = config.tol,
        });

    gelex::notify(dataset_reporter.as_observer(), gelex::DatasetSectionEvent{});

    auto fam_index = gelex::read_fam(config.bfile_prefix + ".fam").index();

    gelex::PhenoPipe pheno(pheno_config, pheno_reporter.as_observer());

    auto grm_paths = std::ranges::to<std::vector<std::filesystem::path>>(
        cmd.get_option("--grm")->as<std::vector<std::string>>());

    gelex::GrmPipe grm(grm_paths, grm_reporter.as_observer());

    auto common = cli::intersect_or_throw(
        dataset_reporter.as_observer(),
        "phenotype, genotype (.fam), GRM, and covariates",
        std::array{&fam_index},
        pheno.sample_indices(),
        grm.sample_indices());

    pheno.load(common);
    grm.load(common);

    cli::RemlReporter reml_reporter;
    gelex::AssocEngine engine(std::move(config));
    engine.run(
        pheno,
        grm,
        common,
        reporter.as_observer(),
        reml_reporter.as_observer());

    return 0;
}
