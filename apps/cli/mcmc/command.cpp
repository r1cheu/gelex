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

#include <cstdint>
#include <filesystem>
#include <optional>
#include <string>
#include <utility>
#include <vector>

#include <Eigen/Core>

#include "cli/bayes_recipe_options.h"
#include "cli/cli_helper.h"
#include "cli/common_data.h"
#include "cli/formatter.h"
#include "cli/geno_reporter.h"
#include "cli/report_printer.h"
#include "gelex/algo/mcmc/params.h"
#include "gelex/algo/mcmc/solver.h"
#include "gelex/bayes/design.h"
#include "gelex/bayes/model.h"
#include "gelex/bayes/recipe.h"
#include "gelex/data/dataframe/index.h"
#include "gelex/data/genotype_method.h"
#include "gelex/data/genotype_reader.h"
#include "gelex/data/reader.h"
#include "gelex/exception.h"
#include "gelex/io/mcmc.h"
#include "gelex/io/snpstats/writer.h"
#include "gelex/types/genetic_mode.h"
#include "reporter.h"

class MCMCDataHandler
{
   public:
    MCMCDataHandler(
        const cli::McmcConfig& config,
        gelex::GeneticModeSet requested_effects,
        cli::GenoReporter& reporter)
        : config_(config),
          requested_effects_(requested_effects),
          reporter_(reporter)
    {
    }

    auto load_indices(std::vector<gelex::DataFrameIndex<std::string>*>& indices)
        -> void
    {
        genotype_method_ = config_.geno_method;

        fam_index_ = gelex::read_fam(config_.bfile + ".fam").index();
        indices.push_back(&fam_index_);
    }

    auto gather(const gelex::DataFrameIndex<std::string>& common_index) -> void
    {
        auto reader = gelex::GenotypeReader(
            config_.bfile, common_index, reporter_.as_observer());

        gelex::BinaryWriter writer(config_.out + ".snpstats");

        for (const auto mode : requested_effects_.each())
        {
            auto genotype
                = config_.mmap
                      ? reader.read_mmap(
                            mode,
                            genotype_method_,
                            std::filesystem::path{
                                config_.out
                                + (mode == gelex::GeneticMode::A ? ".add"
                                                                 : ".dom")},
                            static_cast<std::size_t>(config_.chunk_size))
                      : reader.read_in_memory(
                            mode,
                            genotype_method_,
                            static_cast<std::size_t>(config_.chunk_size));

            reporter_.show_loaded(
                mode, genotype.cols(), genotype.num_invalid());

            gelex::write_snp_stats(writer, mode, genotype.stats());
            genetics_.emplace_back(mode, std::move(genotype));
        }
    }

    auto results() && -> std::vector<gelex::bayes::GeneticDesign>
    {
        return std::move(genetics_);
    }

   private:
    const cli::McmcConfig& config_;
    gelex::GeneticModeSet requested_effects_;
    gelex::GenotypeMethod genotype_method_;
    gelex::DataFrameIndex<std::string> fam_index_;
    std::vector<gelex::bayes::GeneticDesign> genetics_;
    cli::GenoReporter& reporter_;
};

auto mcmc_execute(const cli::McmcConfig& config) -> int
{
    auto recipe_options = cli::make_bayes_recipe_options(config);
    gelex::Params params{
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

    cli::printer().block(gelex::section("[Dataset Summary]"));

    MCMCDataHandler handler(config, recipe_options.modes, geno_reporter);
    cli::BaseData data = cli::load_base_data(handler, config.base_data);
    cli::printer().line(
        "   Intersection : {} common samples", data.sample_ids.size());
    if (data.sample_ids.empty())
    {
        throw gelex::GelexException(
            "No common samples across phenotype, genotype (.fam), and "
            "covariates. Check that sample IDs match across input files.");
    }
    auto bayes_recipe = gelex::bayes::BayesRecipe(std::move(recipe_options));

    auto model = gelex::BayesModel(
        std::move(data.phenotype),
        std::move(data.fixed_design),
        {},
        std::move(handler).results());
    auto prior = bayes_recipe.make_prior(model);

    reporter.show_prior(prior);
    auto result = [&]() -> gelex::Result
    {
        if (config.from_ckpt)
        {
            return solver.run_from(
                model, prior, *config.from_ckpt, reporter.as_observer());
        }
        return solver.run(model, prior, config.seed, reporter.as_observer());
    }();
    reporter.show_complete(result.samples_collected());
    gelex::write_params(result, config.out);
    gelex::write_summary(result, config.out);
    gelex::write_snp_eff(result, model, config.bfile + ".bim", config.out);
    reporter.show_results_saved(config.out);

    return 0;
}
