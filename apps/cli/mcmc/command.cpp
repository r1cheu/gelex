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

#include <Eigen/Core>
#include <cstdint>
#include <optional>
#include <string>
#include <utility>
#include <vector>

#include "gelex/algo/mcmc/params.h"
#include "gelex/algo/mcmc/solver.h"
#include "gelex/bayes/design.h"
#include "gelex/bayes/legacy_recipe.h"
#include "gelex/bayes/model.h"
#include "gelex/data/bed.h"
#include "gelex/data/dataframe/index.h"
#include "gelex/data/genotype_method.h"
#include "gelex/exception.h"
#include "gelex/io/mcmc.h"
#include "gelex/io/snp_lut.h"
#include "gelex/types/genetic_mode.h"

#include "cli/bayes_recipe_options.h"
#include "cli/common_data.h"
#include "cli/formatter.h"
#include "cli/report_printer.h"
#include "cli/runtime.h"
#include "reporter.h"

class MCMCDataHandler
{
   public:
    MCMCDataHandler(cli::GenoReporter& reporter, gelex::Bed& bed)
        : reporter_(reporter), bed_(bed)
    {
    }

    auto load_indices(
        std::vector<const gelex::DataFrameIndex<std::string>*>& indices) -> void
    {
        indices.push_back(&bed_.sample_index());
    }

    auto gather(const gelex::DataFrameIndex<std::string>& common_index) -> void
    {
        bed_.gather(common_index);
        reporter_.show_total(bed_.num_snps());
    }

    auto results() && -> gelex::Bed { return std::move(bed_); }

   private:
    cli::GenoReporter& reporter_;
    gelex::Bed& bed_;
};

auto mcmc_execute(const cli::McmcConfig& config) -> int
{
    auto recipe_options = cli::make_bayes_recipe_options(config);
    gelex::MCMCParams params{
        .n_iters = config.iters,
        .n_burn_in = config.burn_in,
        .n_thin = config.thin,
        .checkpoint_step = config.checkpoint_step.value_or(config.iters),
    };

    gelex::Solver solver{
        params,
        config.out + ".draws",
        params.checkpoint_step > 0 ? std::make_optional(config.out)
                                   : std::nullopt};

    cli::McmcReporter reporter;
    cli::GenoReporter geno_reporter;
    cli::setup_parallelization(config.threads);

    auto bed = gelex::open_bed(config.bfile);
    MCMCDataHandler handler(geno_reporter, bed);

    cli::printer().block(cli::section("Genotype Processing:"));
    cli::BaseData data = cli::load_base_data(handler, config.base_data);
    if (data.sample_ids.empty())
    {
        throw gelex::GelexException(
            "No common samples across phenotype, genotype (.fam), and "
            "covariates. Check that sample IDs match across input files.");
    }

    auto genetic_design = gelex::bayes::GeneticDesign{
        std::move(handler).results(),
        recipe_options.modes,
        config.geno_method,
        geno_reporter.as_observer()};

    auto bayes_recipe = gelex::bayes::BayesRecipe(std::move(recipe_options));

    auto model = gelex::BayesModel(
        std::move(data.phenotype),
        std::move(data.fixed_design),
        {},
        std::move(genetic_design));

    geno_reporter.show_loaded(model.genetic());
    gelex::write_snp_luts(config.out + ".snplut", model.genetic());

    reporter.show_dataset_summary(model, data.pheno_name);
    auto prior = bayes_recipe.make_prior(model);

    reporter.show_prior(prior, model);
    auto result = [&]() -> gelex::Result
    {
        if (config.from_ckpt)
        {
            return solver.run_from(
                model, prior, *config.from_ckpt, reporter.as_observer());
        }
        return solver.run(model, prior, config.seed, reporter.as_observer());
    }();
    reporter.show_summary(result);
    gelex::write_params(result, config.out);
    gelex::write_summary(result, config.out);
    gelex::write_snp_eff(result, model, config.bfile + ".bim", config.out);
    cli::printer().block(
        cli::results_saved(config.out, ".param, .summary, .snpeff, .log"));

    return 0;
}
