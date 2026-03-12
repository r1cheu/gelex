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

#include <optional>

#include <Eigen/Core>

#include <fmt/format.h>

#include "gelex/algo/infer/mcmc.h"
#include "gelex/exception.h"
#include "gelex/infra/logging/notify.h"
#include "gelex/model/bayes/model.h"
#include "gelex/model/bayes/prior_strategies.h"
#include "gelex/model/bayes/trait_model.h"
#include "gelex/model/bayes/writer/result_writer.h"
#include "gelex/pipeline/geno_pipe.h"
#include "gelex/pipeline/pheno_pipe.h"

namespace gelex
{

namespace
{

auto get_default_pi(BayesBase base) -> Eigen::VectorXd
{
    switch (base)
    {
        case BayesBase::B:
        case BayesBase::C:
            return Eigen::VectorXd{{0.99, 0.01}};
        case BayesBase::R:
            return Eigen::VectorXd{{0.99, 0.005, 0.001, 0.001, 0.001}};
        case BayesBase::A:
        case BayesBase::RR:
            return Eigen::VectorXd{{0.0, 1.0}};
    }
    return Eigen::VectorXd{};
}

auto get_default_scale(BayesBase base) -> Eigen::VectorXd
{
    switch (base)
    {
        case BayesBase::R:
            return Eigen::VectorXd{{0.0, 0.001, 0.01, 0.1, 1.0}};
        default:
            return Eigen::VectorXd{};
    }
}

using DefaultFunc = Eigen::VectorXd (*)(BayesBase);

auto to_eigen(
    const std::optional<std::vector<double>>& opt_vec,
    BayesBase base,
    DefaultFunc default_func) -> Eigen::VectorXd
{
    if (opt_vec)
    {
        return Eigen::Map<const Eigen::VectorXd>(
            opt_vec->data(), static_cast<Eigen::Index>(opt_vec->size()));
    }
    return default_func(base);
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

    auto base = config.method.base;
    PriorConfig prior_config;
    prior_config.phenotype_variance = model.phenotype_variance();
    prior_config.positive_prob = config.positive_prob;
    prior_config.genetics.push_back(
        {GeneticEffectType::Add,
         {to_eigen(config.pi, base, get_default_pi),
          to_eigen(config.scale, base, get_default_scale),
          0.5}});
    prior_config.genetics.push_back(
        {GeneticEffectType::Dom,
         {to_eigen(config.dpi, base, get_default_pi),
          to_eigen(config.dscale, base, get_default_scale),
          0.2}});

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

    const auto& m = config.method;
    switch (m.base)
    {
        case BayesBase::A:
            if (m.dominance)
            {
                run_and_write(BayesAd{});
            }
            else
            {
                run_and_write(BayesA{});
            }
            break;
        case BayesBase::B:
            if (m.dominance && m.estimate_pi)
            {
                run_and_write(BayesBdpi{});
            }
            else if (m.dominance)
            {
                run_and_write(BayesBd{});
            }
            else if (m.estimate_pi)
            {
                run_and_write(BayesBpi{});
            }
            else
            {
                run_and_write(BayesB{});
            }
            break;
        case BayesBase::C:
            if (m.dominance && m.estimate_pi)
            {
                run_and_write(BayesCdpi{});
            }
            else if (m.dominance)
            {
                run_and_write(BayesCd{});
            }
            else if (m.estimate_pi)
            {
                run_and_write(BayesCpi{});
            }
            else
            {
                run_and_write(BayesC{});
            }
            break;
        case BayesBase::R:
            if (m.dominance && m.asymmetric)
            {
                run_and_write(BayesRdAt{});
            }
            else if (m.dominance)
            {
                run_and_write(BayesRd{});
            }
            else
            {
                run_and_write(BayesR{});
            }
            break;
        case BayesBase::RR:
            if (m.dominance)
            {
                run_and_write(BayesRRd{});
            }
            else
            {
                run_and_write(BayesRR{});
            }
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
    auto phenotype = std::move(pheno).take_phenotype();
    auto fixed_effects = std::move(pheno).take_fixed_effects();
    auto additive = std::move(geno).take_additive_matrix();
    std::optional<bayes::GenotypeStorage> dominance;
    if (geno.has_dominance_matrix())
    {
        dominance.emplace(std::move(geno).take_dominance_matrix());
    }
    BayesModel model(
        std::move(phenotype),
        std::move(fixed_effects),
        std::move(additive),
        std::move(dominance));
    configure_model_priors(model, config_);

    run_mcmc_analysis(model, config_, observer);

    notify(observer, FitResultsSavedEvent{.out_prefix = config_.out_prefix});
}

}  // namespace gelex
