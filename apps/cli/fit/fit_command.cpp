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

#include "fit_command.h"

#include <argparse.h>
#include <fmt/format.h>
#include <spdlog/logger.h>

#include <Eigen/Core>
#include <filesystem>
#include <memory>
#include <vector>

#include "cli/cli_helper.h"
#include "fit_config.h"
#include "gelex/algo/infer/mcmc.h"
#include "gelex/algo/infer/params.h"
#include "gelex/infra/logger.h"
#include "gelex/infra/utils/formatter.h"
#include "gelex/model/bayes/model.h"
#include "gelex/model/bayes/prior_strategies.h"
#include "gelex/model/bayes/trait_model.h"
#include "gelex/model/effects.h"
#include "gelex/pipeline/data_pipe.h"
#include "gelex/pipeline/report/result_writer.h"

namespace
{

using gelex::BayesAlphabet;

auto get_default_pi(BayesAlphabet type) -> Eigen::VectorXd
{
    switch (type)
    {
        case BayesAlphabet::B:
        case BayesAlphabet::Bpi:
        case BayesAlphabet::Bd:
        case BayesAlphabet::Bdpi:
        case BayesAlphabet::C:
        case BayesAlphabet::Cpi:
        case BayesAlphabet::Cd:
        case BayesAlphabet::Cdpi:
            return Eigen::VectorXd{{0.99, 0.01}};
        case BayesAlphabet::R:
        case BayesAlphabet::Rd:
        {
            return Eigen::VectorXd{{0.99, 0.005, 0.001, 0.001, 0.001}};
        }
        case BayesAlphabet::A:
        case BayesAlphabet::RR:
        case BayesAlphabet::Ad:
        case BayesAlphabet::RRd:
            return Eigen::VectorXd{{0.0, 1.0}};
        default:
            return Eigen::VectorXd{};
    }
}

auto get_default_scale(BayesAlphabet type) -> Eigen::VectorXd
{
    switch (type)
    {
        case BayesAlphabet::R:
        case BayesAlphabet::Rd:
        {
            return Eigen::VectorXd{{0.0, 0.001, 0.01, 0.1, 1.0}};
        }
        default:
            return Eigen::VectorXd{};
    }
}

auto configure_model_priors(
    gelex::BayesModel& model,
    BayesAlphabet type,
    argparse::ArgumentParser& fit,
    const std::shared_ptr<spdlog::logger>& logger) -> int
{
    auto prior_strategy = gelex::create_prior_strategy(type);
    if (!prior_strategy)
    {
        logger->error(
            "Failed to create prior strategy for model type: {}",
            fit.get("-m"));
        return 1;
    }

    gelex::PriorConfig prior_config;
    prior_config.phenotype_variance = model.phenotype_variance();

    auto configure_parameter = [&](std::string_view arg_name,
                                   Eigen::VectorXd& target,
                                   auto default_func)
    {
        if (fit.is_used(arg_name))
        {
            auto values = fit.get<std::vector<double>>(arg_name);
            target = Eigen::Map<const Eigen::VectorXd>(
                values.data(), static_cast<Eigen::Index>(values.size()));
        }
        else
        {
            target = default_func(type);
        }
    };

    configure_parameter(
        "--pi", prior_config.additive.mixture_proportions, get_default_pi);
    configure_parameter(
        "--dpi", prior_config.dominant.mixture_proportions, get_default_pi);
    configure_parameter(
        "--scale", prior_config.additive.mixture_scales, get_default_scale);
    configure_parameter(
        "--dscale", prior_config.dominant.mixture_scales, get_default_scale);

    (*prior_strategy)(model, prior_config);
    return 0;
}

auto run_mcmc_analysis(
    gelex::BayesModel& model,
    BayesAlphabet type,
    int seed,
    const gelex::MCMCParams& mcmc_params,
    const std::filesystem::path& bim_path,
    std::string_view out_prefix,
    const std::shared_ptr<spdlog::logger>& logger) -> int
{
    auto run_and_write = [&](auto trait_model)
    {
        gelex::MCMC mcmc(mcmc_params, trait_model);
        gelex::MCMCResult result = mcmc.run(model, seed, out_prefix);
        gelex::MCMCResultWriter writer(result, bim_path);
        writer.save(out_prefix);
    };

    switch (type)
    {
        case (BayesAlphabet::A):
            run_and_write(gelex::BayesA{});
            break;
        case (BayesAlphabet::Ad):
            run_and_write(gelex::BayesAd{});
            break;
        case (BayesAlphabet::B):
            run_and_write(gelex::BayesB{});
            break;
        case (BayesAlphabet::Bpi):
            run_and_write(gelex::BayesBpi{});
            break;
        case (BayesAlphabet::Bd):
            run_and_write(gelex::BayesBd{});
            break;
        case (BayesAlphabet::Bdpi):
            run_and_write(gelex::BayesBdpi{});
            break;
        case (BayesAlphabet::C):
            run_and_write(gelex::BayesC{});
            break;
        case (BayesAlphabet::Cpi):
            run_and_write(gelex::BayesCpi{});
            break;
        case (BayesAlphabet::Cd):
            run_and_write(gelex::BayesCd{});
            break;
        case (BayesAlphabet::Cdpi):
            run_and_write(gelex::BayesCdpi{});
            break;
        case (BayesAlphabet::R):
            run_and_write(gelex::BayesR{});
            break;
        case (BayesAlphabet::Rd):
            run_and_write(gelex::BayesRd{});
            break;
        case (BayesAlphabet::RR):
            run_and_write(gelex::BayesRR{});
            break;
        case (BayesAlphabet::RRd):
            run_and_write(gelex::BayesRRd{});
            break;
        default:
            logger->error("Unsupported method: {}", type);
            return 1;
    }
    return 0;
}

}  // namespace

int fit_execute(argparse::ArgumentParser& fit)
{
    auto config = FitConfig::make(fit);

    auto logger = gelex::logging::get();

    gelex::cli::setup_parallelization(config.threads);

    gelex::cli::print_fit_header(
        config.method_name,
        config.use_dominance,
        static_cast<int>(config.mcmc_params.n_iters),
        static_cast<int>(config.mcmc_params.n_burnin),
        config.threads);

    // ================================================================
    // Data Loading & Pipeline
    // ================================================================
    gelex::DataPipe::Config data_pipe_config{
        .phenotype_path = config.phenotype_path,
        .phenotype_column = config.phenotype_column,
        .bed_path = config.bed_path,
        .use_dominance_effect = config.use_dominance,
        .use_mmap = config.use_mmap,
        .chunk_size = config.chunk_size,
        .qcovar_path = config.qcovar_path,
        .dcovar_path = config.dcovar_path,
        .output_prefix = config.out_prefix,
        .genotype_method = config.genotype_method};

    gelex::DataPipe data_pipe(data_pipe_config);
    logger->info(gelex::section("[Dataset Summary]"));
    auto p_stats = data_pipe.load_phenotypes();
    logger->info(
        gelex::success(
            "Phenotypes : {} samples ('{}')",
            p_stats.samples_loaded,
            p_stats.trait_name));
    logger->info(
        gelex::success(
            "Genotypes  : {} samples", data_pipe.num_genotype_samples()));

    auto c_stats = data_pipe.load_covariates();
    if (c_stats.qcovar_loaded > 0 || c_stats.dcovar_loaded > 0)
    {
        std::string parts;
        if (c_stats.qcovar_loaded > 0)
        {
            parts += fmt::format(
                "{} quantitative ({})",
                c_stats.qcovar_loaded,
                gelex::format_names(c_stats.q_names));
        }
        if (c_stats.qcovar_loaded > 0 && c_stats.dcovar_loaded > 0)
        {
            parts += ", ";
        }
        if (c_stats.dcovar_loaded > 0)
        {
            parts += fmt::format(
                "{} discrete ({})",
                c_stats.dcovar_loaded,
                gelex::format_names(c_stats.d_names));
        }
        logger->info(gelex::success("Covariates : {}", parts));
    }
    auto i_stats = data_pipe.intersect_samples();
    logger->info(
        gelex::success(
            "Intersection : {} common, {} excluded",
            i_stats.common_samples,
            i_stats.excluded_samples));

    if (i_stats.common_samples == 0)
    {
        logger->error(
            "No common samples found between phenotype, covariates, and "
            "genotype files.");
        return 1;
    }

    auto log_overwrite = [&](const std::string& msg)
    {
        if (gelex::cli::is_tty())
        {
            logger->info("{}", "\033[A\r" + msg + "\033[K");
        }
        else
        {
            logger->info("{}", msg);
        }
    };

    auto add_stats = data_pipe.load_additive_matrix();
    log_overwrite(
        gelex::success(
            "Additive     : {} SNPs ({} monomorphic excluded)",
            gelex::AbbrNumber(add_stats.num_snps - add_stats.monomorphic_snps),
            gelex::AbbrNumber(add_stats.monomorphic_snps)));
    if (config.use_dominance)
    {
        auto dom_stats = data_pipe.load_dominance_matrix();
        log_overwrite(
            gelex::success(
                "Dominance    : {} SNPs ({} monomorphic excluded)",
                gelex::AbbrNumber(
                    add_stats.num_snps - add_stats.monomorphic_snps),
                gelex::AbbrNumber(dom_stats.monomorphic_snps)));
    }

    data_pipe.finalize();

    gelex::BayesModel model(data_pipe);

    if (configure_model_priors(model, config.method, fit, logger) != 0)
    {
        return 1;
    }

    if (run_mcmc_analysis(
            model,
            config.method,
            config.seed,
            config.mcmc_params,
            config.bed_path.replace_extension(".bim"),
            config.out_prefix,
            logger)
        != 0)
    {
        return 1;
    }
    logger->info(
        gelex::success(
            "Results saved to '{}' (.param, .snp.eff, .log)",
            config.out_prefix));

    return 0;
}
