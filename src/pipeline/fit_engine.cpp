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

#include <Eigen/Core>

#include <fmt/format.h>

#include "gelex/algo/infer/mcmc.h"
#include "gelex/exception.h"
#include "gelex/infra/logging/notify.h"
#include "gelex/model/bayes/model.h"
#include "gelex/model/bayes/prior_strategies.h"
#include "gelex/model/bayes/trait_model.h"
#include "gelex/pipeline/geno_pipe.h"
#include "gelex/pipeline/pheno_pipe.h"
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

auto configure_model_priors(BayesModel& model, const FitEngine::Config& config)
    -> void
{
    auto prior_strategy = create_prior_strategy(config.method);
    if (!prior_strategy)
    {
        throw GelexException(
            fmt::format(
                "Failed to create prior strategy for model type: {}",
                config.method));
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
        auto bim_path = config.bfile_prefix + ".bim";
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

auto FitEngine::run(
    PhenoPipe&& pheno,
    GenoPipe&& geno,
    const FitObserver& observer) -> void
{
    auto pheno_pipe = std::move(pheno);
    auto geno_pipe = std::move(geno);
    BayesModel model(pheno_pipe, geno_pipe);
    configure_model_priors(model, config_);

    run_mcmc_analysis(model, config_, observer);

    notify(observer, FitResultsSavedEvent{.out_prefix = config_.out_prefix});
}

}  // namespace gelex
