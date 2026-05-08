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

#include "assoc_command.h"

#include <argparse.h>
#include <filesystem>
#include <ranges>
#include <string>
#include <vector>

#include "assoc_config.h"
#include "assoc_reporter.h"
#include "cli/cli_helper.h"
#include "cli/data_pipe_config.h"
#include "cli/dataset_loader.h"
#include "cli/dataset_reporter.h"
#include "cli/grm_pipe_reporter.h"
#include "cli/pheno_reporter.h"
#include "cli/reml_reporter.h"
#include "gelex/data/dataframe/index.h"
#include "gelex/data/genotype/bed_path.h"
#include "gelex/data/pipe/grm.h"
#include "gelex/data/pipe/pheno.h"
#include "gelex/data/reader.h"
#include "gelex/engine/assoc.h"
#include "gelex/infra/logging/assoc_event.h"
#include "gelex/infra/logging/dataset_event.h"
#include "gelex/infra/logging/notify.h"

auto assoc_execute(argparse::ArgumentParser& cmd) -> int
{
    int threads = cmd.get<int>("--threads");
    gelex::cli::setup_parallelization(threads);

    auto pheno_config = gelex::cli::make_pheno_config(cmd);
    pheno_config.transform_type
        = gelex::cli::parse_transform_type(cmd.get("--transform"));
    pheno_config.int_offset = cmd.get<double>("--int-offset");

    gelex::cli::AssocReporter reporter;
    gelex::cli::DatasetReporter dataset_reporter;
    gelex::cli::PhenoReporter pheno_reporter;
    gelex::cli::GrmPipeReporter grm_reporter;

    auto config = gelex::cli::make_assoc_config(cmd);

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

    auto bed_path
        = gelex::genotype::format_bed_path(cmd.get<std::string>("--bfile"));
    auto fam_index
        = gelex::read_fam(
              std::filesystem::path(bed_path).replace_extension(".fam"))
              .index();

    gelex::PhenoPipe pheno(pheno_config, pheno_reporter.as_observer());

    auto grm_paths = std::ranges::to<std::vector<std::filesystem::path>>(
        cmd.get<std::vector<std::string>>("--grm"));

    gelex::GrmPipe grm(grm_paths, grm_reporter.as_observer());

    std::vector<const gelex::dataframe::Index<std::string>*> all_indices{
        &fam_index};
    all_indices.append_range(pheno.sample_indices());
    all_indices.append_range(grm.sample_indices());
    auto common = gelex::cli::intersect_or_throw(
        std::move(all_indices),
        dataset_reporter.as_observer(),
        "phenotype, genotype (.fam), GRM, and covariates");

    pheno.load(common);
    grm.load(common);

    gelex::cli::RemlReporter reml_reporter;
    gelex::AssocEngine engine(std::move(config));
    engine.run(
        pheno,
        grm,
        common,
        reporter.as_observer(),
        reml_reporter.as_observer());

    return 0;
}
