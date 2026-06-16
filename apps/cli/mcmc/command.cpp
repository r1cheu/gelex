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

#include <CLI/CLI.hpp>
#include <Eigen/Core>

#include "cli/bayes_recipe_options.h"
#include "cli/cli_helper.h"
#include "cli/common_data.h"
#include "cli/dataset_reporter.h"
#include "cli/geno_reporter.h"
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
#include "gelex/infra/logging/dataset_event.h"
#include "gelex/infra/logging/fit_event.h"
#include "gelex/infra/logging/geno_event.h"
#include "gelex/infra/logging/notify.h"
#include "gelex/io/locistats/writer.h"
#include "gelex/io/mcmc.h"
#include "gelex/types/genetic_effect_type.h"
#include "reporter.h"

class MCMCDataHandler
{
   public:
    MCMCDataHandler(
        std::vector<gelex::GeneticMode> requested_effects,
        gelex::GenoObserver observer)
        : requested_effects_(std::move(requested_effects)),
          observer_(std::move(observer))
    {
    }

    auto load_indices(
        CLI::App& cmd,
        std::vector<gelex::dataframe::Index<std::string>*>& indices) -> void
    {
        bfile_prefix_ = cmd.get_option("--bfile")->as<std::string>();
        out_prefix_ = cmd.get_option("--out")->as<std::string>();
        chunk_size_ = cmd.get_option("--chunk-size")->as<int>();
        if (chunk_size_ <= 0)
        {
            throw gelex::GelexException("--chunk-size must be positive");
        }
        genotype_method_ = cli::parse_genotype_method(
            cmd.get_option("--geno-method")->as<std::string>());
        use_mmap_ = cmd.get_option("--mmap")->count() > 0;

        fam_index_ = gelex::read_fam(bfile_prefix_ + ".fam").index();
        indices.push_back(&fam_index_);
    }

    auto gather(const gelex::dataframe::Index<std::string>& common_index)
        -> void
    {
        auto reader = gelex::genotype::GenotypeReader(
            bfile_prefix_, common_index, observer_);

        gelex::LociStatsWriter writer(out_prefix_ + ".sbin");
        const auto method_code
            = static_cast<std::uint8_t>(std::to_underlying(genotype_method_));
        const bool method_is_center = gelex::is_center(genotype_method_);

        for (const auto mode : requested_effects_)
        {
            auto genotype
                = use_mmap_
                      ? reader.read_mmap(
                            mode,
                            genotype_method_,
                            std::filesystem::path{
                                out_prefix_
                                + (mode == gelex::GeneticMode::A ? ".add"
                                                                 : ".dom")},
                            static_cast<std::size_t>(chunk_size_))
                      : reader.read_in_memory(
                            mode,
                            genotype_method_,
                            static_cast<std::size_t>(chunk_size_));

            gelex::notify(
                observer_,
                gelex::GenotypeLoadedEvent{
                    .mode = mode,
                    .num_snps = genotype.cols(),
                    .monomorphic_snps = genotype.num_mono()});

            const Eigen::VectorXd stddev{genotype.var().array().sqrt()};
            writer.write(
                gelex::EffectType::from_genetic(mode),
                method_code,
                genotype.mean(),
                method_is_center ? static_cast<const Eigen::VectorXd*>(nullptr)
                                 : &stddev,
                genotype.mono_indices());
            genetics_.emplace_back(mode, std::move(genotype));
        }
    }

    auto results() && -> std::vector<gelex::bayes::GeneticDesign>
    {
        return std::move(genetics_);
    }

   private:
    std::vector<gelex::GeneticMode> requested_effects_;
    std::string bfile_prefix_;
    std::string out_prefix_;
    gelex::GenotypeMethod genotype_method_;
    int chunk_size_;
    bool use_mmap_;
    gelex::dataframe::Index<std::string> fam_index_;
    std::vector<gelex::bayes::GeneticDesign> genetics_;
    gelex::GenoObserver observer_;
};

auto mcmc_execute(CLI::App& cmd) -> int
{
    auto recipe_options = cli::make_bayes_recipe_options(cmd);
    const auto n_iters = cmd.get_option("--iters")->as<int>();
    gelex::mcmc::Params params{
        .n_iters = n_iters,
        .n_burn_in = cmd.get_option("--burn-in")->as<int>(),
        .n_thin = cmd.get_option("--thin")->as<int>(),
        .checkpoint_step = cmd.get_option("--checkpoint-step")->count() > 0
                               ? cmd.get_option("--checkpoint-step")->as<int>()
                               : n_iters,
    };

    gelex::mcmc::Solver solver{
        params,
        cmd.get_option("--out")->as<std::string>() + ".draws",
        params.checkpoint_step > 0
            ? std::make_optional(cmd.get_option("--out")->as<std::string>())
            : std::nullopt};

    const int threads = cmd.get_option("--threads")->as<int>();
    cli::McmcReporter reporter;
    cli::DatasetReporter dataset_reporter;
    cli::GenoReporter geno_reporter;
    cli::setup_parallelization(threads);

    gelex::notify(reporter.as_observer(), gelex::MCMCBannerEvent{});
    gelex::notify(
        reporter.as_observer(),
        gelex::MCMCConfigEvent{
            .recipe_scheme = recipe_options.scheme,
            .requested_effects = recipe_options.modes,
            .n_iters = params.n_iters,
            .n_burn_in = params.n_burn_in,
            .seed = cmd.get_option("--seed")->as<int>(),
        });

    gelex::notify(dataset_reporter.as_observer(), gelex::DatasetSectionEvent{});

    MCMCDataHandler handler(recipe_options.modes, geno_reporter.as_observer());
    cli::BaseData data = cli::load_base_data(handler, cmd);
    gelex::notify(
        dataset_reporter.as_observer(),
        gelex::IntersectionEvent{.common_samples = data.sample_ids.size()});
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

    gelex::notify(reporter.as_observer(), gelex::FitPriorSetEvent{&prior});
    auto result = [&]() -> gelex::mcmc::Result
    {
        if (cmd.get_option("--from-ckpt")->count() > 0)
        {
            return solver.run_from(
                model,
                prior,
                cmd.get_option("--from-ckpt")->as<std::string>(),
                reporter.as_observer());
        }
        return solver.run(
            model,
            prior,
            cmd.get_option("--seed")->as<int>(),
            reporter.as_observer());
    }();
    gelex::mcmc::write_params(
        result, cmd.get_option("--out")->as<std::string>());
    gelex::mcmc::write_summary(
        result, cmd.get_option("--out")->as<std::string>());
    gelex::mcmc::write_snp_eff(
        result,
        model,
        cmd.get_option("--bfile")->as<std::string>() + ".bim",
        cmd.get_option("--out")->as<std::string>());
    gelex::notify(
        reporter.as_observer(),
        gelex::FitResultsSavedEvent{
            .out_prefix = cmd.get_option("--out")->as<std::string>()});

    return 0;
}
