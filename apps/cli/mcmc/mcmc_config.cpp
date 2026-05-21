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

#include <array>
#include <optional>
#include <string_view>
#include <vector>

#include <argparse.h>
#include <Eigen/Core>

#include <fmt/format.h>

#include "cli/cli_helper.h"
#include "gelex/exception.h"
#include "gelex/types/genetic_effect_type.h"

namespace gelex::cli
{

namespace
{

auto reject_canonical_bayes_recipe_flags(const argparse::ArgumentParser& cmd)
    -> void
{
    constexpr std::array kFlags = {
        std::string_view{"--additive-h2"},
        std::string_view{"--dominance-h2"},
        std::string_view{"--additive-pi"},
        std::string_view{"--dominance-pi"},
        std::string_view{"--additive-multiplier"},
        std::string_view{"--dominance-multiplier"},
        std::string_view{"--joint-pi"},
        std::string_view{"--estimate-additive-pi"},
        std::string_view{"--estimate-dominance-pi"},
        std::string_view{"--estimate-joint-pi"},
        std::string_view{"--dominance-positive-prob"},
        std::string_view{"--random-variance-proportion"},
    };
    for (const auto flag : kFlags)
    {
        if (cmd.is_used(flag))
        {
            throw GelexException(
                fmt::format(
                    "BayesRecipe flag {} is not supported by the current "
                    "MCMC engine; remove it or use a legacy equivalent where "
                    "available",
                    flag));
        }
    }
}

auto parse_method(argparse::ArgumentParser& cmd)
    -> std::pair<bayes::LegacyBayesConfig, std::vector<GeneticMode>>
{
    auto base
        = gelex::get_bayes_base(cmd.get("-m")).value_or(gelex::BayesBase::RR);
    auto requested = parse_genetic_modes(cmd.get("--mode"));

    bayes::LegacyBayesConfig cfg{
        .base = base,
        .dominance = cmd.get<bool>("--asym")
                         ? bayes::DominancePolicy::asymmetric
                         : bayes::DominancePolicy::symmetric,
        .estimate_pi = cmd.get<bool>("--estimate-pi"),
    };
    return {cfg, std::move(requested)};
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
    bayes::LegacyBayesConfig method,
    std::vector<GeneticMode> requested) -> mcmc::Engine::Config
{
    mcmc::Engine::Config config{
        .bfile_prefix = cmd.get("--bfile"),
        .method = method,
        .requested_effects = std::move(requested),
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

auto make_mcmc_config(argparse::ArgumentParser& cmd) -> McmcConfig
{
    reject_canonical_bayes_recipe_flags(cmd);
    auto [method, requested] = parse_method(cmd);
    auto overrides = make_overrides(cmd);
    validate_overrides(method, overrides);
    return McmcConfig{
        .engine = make_engine_config(cmd, method, std::move(requested)),
        .overrides = std::move(overrides),
    };
}

}  // namespace gelex::cli
