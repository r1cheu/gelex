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

#include <argparse.h>

#include "gelex/exception.h"
#include "gelex/pipeline/fit_engine.h"

namespace gelex::cli
{

auto make_fit_config(argparse::ArgumentParser& cmd) -> FitEngine::Config
{
    auto base
        = gelex::get_bayes_base(cmd.get("-m")).value_or(gelex::BayesBase::RR);

    auto method = gelex::BayesMethodConfig{
        .base = base,
        .dominance = cmd.get<bool>("--dom"),
        .asymmetric = cmd.get<bool>("--asym"),
        .estimate_pi = cmd.get<bool>("--estimate-pi"),
    };

    if (!gelex::is_valid_method(method))
    {
        throw gelex::InvalidInputException(
            fmt::format("invalid method combination: {}", method));
    }

    FitEngine::Config config{
        .bfile_prefix = cmd.get("--bfile"),
        .method = method,

        .seed = cmd.get<int>("--seed"),
        .mcmc_params = gelex::MCMCParams(
            cmd.get<int>("--iters"),
            cmd.get<int>("--burn-in"),
            cmd.get<int>("--thin"),
            cmd.is_used("--checkpoint-step") ? cmd.get<int>("--checkpoint-step")
                                             : 0),
        .out_prefix = cmd.get("--out")};

    auto extract_opt_vec
        = [&](std::string_view arg) -> std::optional<std::vector<double>>
    {
        if (cmd.is_used(arg))
        {
            return cmd.get<std::vector<double>>(arg);
        }
        return std::nullopt;
    };

    config.pi = extract_opt_vec("--pi");
    config.dpi = extract_opt_vec("--dpi");
    config.multiplier = extract_opt_vec("--mult");
    config.dmultiplier = extract_opt_vec("--dmult");
    config.positive_prob = cmd.get<double>("--positive-prob");
    if (config.mcmc_params.n_burn_in >= config.mcmc_params.n_iters)
    {
        throw gelex::InvalidInputException(
            "n_burn_in must be less than n_iters");
    }

    if (cmd.is_used("--resume"))
    {
        config.resume_path = cmd.get("--resume");
    }

    return config;
}

}  // namespace gelex::cli
