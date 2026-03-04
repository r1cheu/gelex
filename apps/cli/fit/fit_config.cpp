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

auto has_dominance(BayesAlphabet type) -> bool
{
    switch (type)
    {
        case BayesAlphabet::Bd:
        case BayesAlphabet::Bdpi:
        case BayesAlphabet::Cd:
        case BayesAlphabet::Cdpi:
        case BayesAlphabet::Rd:
        case BayesAlphabet::Ad:
        case BayesAlphabet::RRd:
            return true;
        default:
            return false;
    }
}

auto make_fit_config(argparse::ArgumentParser& cmd) -> FitEngine::Config
{
    auto method = gelex::get_bayesalphabet(cmd.get("-m"))
                      .value_or(gelex::BayesAlphabet::RR);
    auto use_dominance = has_dominance(method);

    FitEngine::Config config{
        .bfile_prefix = cmd.get("--bfile"),
        .method = method,

        .seed = cmd.get<int>("--seed"),
        .mcmc_params = gelex::MCMCParams(
            cmd.get<int>("--iters"),
            cmd.get<int>("--burn-in"),
            cmd.get<int>("--thin")),
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
    config.scale = extract_opt_vec("--scale");
    config.dscale = extract_opt_vec("--dscale");
    if (config.mcmc_params.n_burnin >= config.mcmc_params.n_iters)
    {
        throw gelex::InvalidInputException(
            "n_burnin must be less than n_iters");
    }

    return config;
}

}  // namespace gelex::cli
