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

#include "mcmc_config.h"

#include <optional>
#include <vector>

#include <argparse.h>
#include <Eigen/Core>

#include "gelex/exception.h"

namespace gelex::cli
{

namespace
{

auto parse_method(argparse::ArgumentParser& cmd) -> bayes::BayesConfig
{
    auto base
        = gelex::get_bayes_base(cmd.get("-m")).value_or(gelex::BayesBase::RR);

    auto mode = gelex::get_genetic_mode(cmd.get("--mode"));
    if (!mode)
    {
        throw gelex::GelexException(
            fmt::format("invalid --mode: {}", cmd.get("--mode")));
    }

    auto method = bayes::BayesConfig{
        .base = base,
        .mode = *mode,
        .dominance = cmd.get<bool>("--asym")
                         ? bayes::DominancePolicy::asymmetric
                         : bayes::DominancePolicy::symmetric,
        .estimate_pi = cmd.get<bool>("--estimate-pi"),
    };

    if (!bayes::is_valid_method(method))
    {
        throw gelex::GelexException(
            fmt::format("invalid method combination: {}", method));
    }
    return method;
}

auto extract_eigen(argparse::ArgumentParser& cmd, std::string_view arg)
    -> std::optional<Eigen::VectorXd>
{
    if (!cmd.is_used(arg))
    {
        return std::nullopt;
    }
    auto v = cmd.get<std::vector<double>>(arg);
    return Eigen::Map<const Eigen::VectorXd>(
        v.data(), static_cast<Eigen::Index>(v.size()));
}

auto make_overrides(argparse::ArgumentParser& cmd) -> MethodOverrides
{
    MethodOverrides o;
    o.additive.proportions = extract_eigen(cmd, "--pi");
    o.additive.multiplier = extract_eigen(cmd, "--mult");
    o.dominance.proportions = extract_eigen(cmd, "--dpi");
    o.dominance.multiplier = extract_eigen(cmd, "--dmult");
    if (cmd.is_used("--positive-prob"))
    {
        o.dominance.positive_prob = cmd.get<double>("--positive-prob");
    }
    return o;
}

auto make_engine_config(
    argparse::ArgumentParser& cmd,
    bayes::BayesConfig method) -> mcmc::FitEngine::Config
{
    mcmc::FitEngine::Config config{
        .bfile_prefix = cmd.get("--bfile"),
        .method = method,
        .seed = cmd.get<int>("--seed"),
        .mcmc_params = gelex::mcmc::Params{
            .n_iters = cmd.get<int>("--iters"),
            .n_burn_in = cmd.get<int>("--burn-in"),
            .n_thin = cmd.get<int>("--thin"),
            .checkpoint_step
            = cmd.is_used("--checkpoint-step")
                  ? cmd.get<int>("--checkpoint-step")
                  : 0,
        },
        .out_prefix = cmd.get("--out"),
    };

    if (cmd.is_used("--resume"))
    {
        config.resume_path = cmd.get("--resume");
    }

    mcmc::ConfigValidator{config}.validate();
    return config;
}

}  // namespace

auto make_mcmc_fit_config(argparse::ArgumentParser& cmd) -> McmcFitConfig
{
    auto method = parse_method(cmd);
    auto overrides = make_overrides(cmd);
    validate_overrides(method, overrides);
    return McmcFitConfig{
        .engine = make_engine_config(cmd, method),
        .overrides = std::move(overrides),
    };
}

}  // namespace gelex::cli
