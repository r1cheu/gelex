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

#include <string>
#include <utility>

#include <argparse.h>

#include "cli/bayes_recipe_config.h"
#include "gelex/algo/infer/params.h"
#include "gelex/bayes/recipe.h"
#include "gelex/bayes/recipe_options.h"
#include "gelex/engine/mcmc.h"

namespace gelex::cli
{

namespace
{

auto make_engine_config(
    argparse::ArgumentParser& cmd,
    bayes::BayesRecipeScheme recipe_scheme,
    bayes::BayesRecipeConfig recipe_config) -> mcmc::Engine::Config
{
    mcmc::Engine::Config config{
        .bfile_prefix = cmd.get("--bfile"),
        .recipe_scheme = recipe_scheme,
        .recipe_config = std::move(recipe_config),
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

auto make_mcmc_engine_config(argparse::ArgumentParser& cmd)
    -> mcmc::Engine::Config
{
    auto recipe_scheme = gelex::bayes::to_bayes_recipe_scheme(
        cmd.get<std::string>("--method"));
    auto recipe_config = make_bayes_recipe_config(cmd);
    return make_engine_config(cmd, recipe_scheme, std::move(recipe_config));
}

}  // namespace gelex::cli
