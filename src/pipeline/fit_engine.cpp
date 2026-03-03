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

#include "gelex/pipeline/fit_engine.h"

#include <unistd.h>

#include <fmt/format.h>
#include <Eigen/Core>

#include <memory>

#include "gelex/algo/infer/mcmc.h"
#include "gelex/infra/logger.h"
#include "gelex/infra/utils/formatter.h"
#include "gelex/model/bayes/model.h"
#include "gelex/model/bayes/prior_strategies.h"
#include "gelex/model/bayes/trait_model.h"
#include "gelex/pipeline/data_pipe.h"
#include "gelex/pipeline/report/result_writer.h"

namespace gelex
{

namespace
{

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
            return Eigen::VectorXd{{0.99, 0.005, 0.001, 0.001, 0.001}};
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
            return Eigen::VectorXd{{0.0, 0.001, 0.01, 0.1, 1.0}};
        default:
            return Eigen::VectorXd{};
    }
}

auto to_eigen(
    const std::optional<std::vector<double>>& opt_vec,
    BayesAlphabet type,
    Eigen::VectorXd (*default_func)(BayesAlphabet)) -> Eigen::VectorXd
{
    if (opt_vec)
    {
        return Eigen::Map<const Eigen::VectorXd>(
            opt_vec->data(), static_cast<Eigen::Index>(opt_vec->size()));
    }
    return default_func(type);
}

auto configure_model_priors(
    BayesModel& model,
    const FitEngine::Config& config,
    const std::shared_ptr<spdlog::logger>& logger) -> void
{
    auto prior_strategy = create_prior_strategy(config.method);
    if (!prior_strategy)
    {
        logger->error(
            "Failed to create prior strategy for model type: {}",
            config.method);
        return;
    }

    PriorConfig prior_config;
    prior_config.phenotype_variance = model.phenotype_variance();
    prior_config.additive.mixture_proportions
        = to_eigen(config.pi, config.method, get_default_pi);
    prior_config.dominant.mixture_proportions
        = to_eigen(config.dpi, config.method, get_default_pi);
    prior_config.additive.mixture_scales
        = to_eigen(config.scale, config.method, get_default_scale);
    prior_config.dominant.mixture_scales
        = to_eigen(config.dscale, config.method, get_default_scale);

    (*prior_strategy)(model, prior_config);
}

auto notify(const FitObserver& observer, const FitEvent& event) -> void
{
    if (observer)
    {
        observer(event);
    }
}

auto run_mcmc_analysis(
    BayesModel& model,
    const FitEngine::Config& config,
    const FitObserver& observer) -> void
{
    auto run_and_write = [&](auto trait_model)
    {
        MCMC mcmc(config.mcmc_params, trait_model);
        MCMCResult result
            = mcmc.run(model, config.seed, config.out_prefix, observer);
        auto bim_path = config.bed_path;
        bim_path.replace_extension(".bim");
        MCMCResultWriter writer(result, bim_path);
        writer.save(config.out_prefix);
    };

    switch (config.method)
    {
        case BayesAlphabet::A:
            run_and_write(BayesA{});
            break;
        case BayesAlphabet::Ad:
            run_and_write(BayesAd{});
            break;
        case BayesAlphabet::B:
            run_and_write(BayesB{});
            break;
        case BayesAlphabet::Bpi:
            run_and_write(BayesBpi{});
            break;
        case BayesAlphabet::Bd:
            run_and_write(BayesBd{});
            break;
        case BayesAlphabet::Bdpi:
            run_and_write(BayesBdpi{});
            break;
        case BayesAlphabet::C:
            run_and_write(BayesC{});
            break;
        case BayesAlphabet::Cpi:
            run_and_write(BayesCpi{});
            break;
        case BayesAlphabet::Cd:
            run_and_write(BayesCd{});
            break;
        case BayesAlphabet::Cdpi:
            run_and_write(BayesCdpi{});
            break;
        case BayesAlphabet::R:
            run_and_write(BayesR{});
            break;
        case BayesAlphabet::Rd:
            run_and_write(BayesRd{});
            break;
        case BayesAlphabet::RR:
            run_and_write(BayesRR{});
            break;
        case BayesAlphabet::RRd:
            run_and_write(BayesRRd{});
            break;
        default:
            break;
    }
}

}  // namespace

FitEngine::FitEngine(Config config) : config_(std::move(config)) {}

auto FitEngine::run(const FitObserver& observer) -> void
{
    auto logger = gelex::logging::get();

    DataPipe::Config data_pipe_config{
        .phenotype_path = config_.phenotype_path,
        .phenotype_column = config_.phenotype_column,
        .bed_path = config_.bed_path,
        .use_dominance_effect = config_.use_dominance,
        .use_mmap = config_.use_mmap,
        .chunk_size = config_.chunk_size,
        .qcovar_path = config_.qcovar_path,
        .dcovar_path = config_.dcovar_path,
        .output_prefix = config_.out_prefix,
        .grm_paths = {},
        .genotype_method = config_.genotype_method,
    };

    DataPipe data_pipe(data_pipe_config);

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
        return;
    }

    auto log_overwrite = [&](const std::string& msg)
    {
        if (isatty(fileno(stdout)) != 0)
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

    if (config_.use_dominance)
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

    BayesModel model(data_pipe);
    configure_model_priors(model, config_, logger);

    run_mcmc_analysis(model, config_, observer);

    notify(observer, FitResultsSavedEvent{.out_prefix = config_.out_prefix});
}

}  // namespace gelex
