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

#include "fit_config.h"

#include <initializer_list>
#include <optional>
#include <string_view>
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

auto reject_if_used(
    argparse::ArgumentParser& cmd,
    std::initializer_list<std::string_view> flags,
    std::string_view context) -> void
{
    for (auto flag : flags)
    {
        if (cmd.is_used(flag))
        {
            throw gelex::GelexException(
                fmt::format("{} is not valid with {}", flag, context));
        }
    }
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

auto make_mcmc_config(argparse::ArgumentParser& cmd, bayes::BayesConfig method)
    -> mcmc::FitEngine::Config
{
    reject_if_used(cmd, {"--max-iters", "--tol"}, "--im mcmc");

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

auto make_cavi_config(argparse::ArgumentParser& cmd, bayes::BayesConfig method)
    -> vi::FitEngine::Config
{
    if (method.base != BayesBase::RR)
    {
        throw gelex::GelexException(
            fmt::format("CAVI only supports BayesRR, got: {}", method));
    }

    reject_if_used(
        cmd,
        {"--iters",
         "--burn-in",
         "--thin",
         "--seed",
         "--checkpoint-step",
         "--resume",
         "--pi",
         "--dpi",
         "--mult",
         "--dmult",
         "--estimate-pi",
         "--asym",
         "--positive-prob"},
        "--im cavi");

    return vi::FitEngine::Config{
        .bfile_prefix = cmd.get("--bfile"),
        .method = method,
        .params = vi::Params{
            .max_iters = cmd.get<int>("--max-iters"),
            .tol = cmd.get<double>("--tol"),
        },
        .out_prefix = cmd.get("--out"),
    };
}

}  // namespace

auto make_fit_config(argparse::ArgumentParser& cmd) -> FitConfig
{
    auto method = parse_method(cmd);
    auto infer = cmd.get("--infer-method");

    if (infer == "cavi")
    {
        return FitConfig{.engine = make_cavi_config(cmd, method)};
    }

    auto overrides = make_overrides(cmd);
    validate_overrides(method, overrides);
    return FitConfig{
        .engine = make_mcmc_config(cmd, method),
        .overrides = std::move(overrides),
    };
}

}  // namespace gelex::cli
