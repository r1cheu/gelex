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

#include "config.h"

#include <string>

#include <CLI/CLI.hpp>

#include "gelex/algo/mcmc/params.h"
#include "gelex/engine/mcmc.h"

namespace cli
{

auto make_mcmc_engine_config(CLI::App& cmd) -> gelex::mcmc::Engine::Config
{
    const auto n_iters = cmd.get_option("--iters")->as<int>();
    gelex::mcmc::Engine::Config config{
        .bfile_prefix = cmd.get_option("--bfile")->as<std::string>(),
        .seed = cmd.get_option("--seed")->as<int>(),
        .mcmc_params = gelex::mcmc::Params{
            .n_iters = n_iters,
            .n_burn_in = cmd.get_option("--burn-in")->as<int>(),
            .n_thin = cmd.get_option("--thin")->as<int>(),
            .checkpoint_step
            = cmd.get_option("--checkpoint-step")->count() > 0
                  ? cmd.get_option("--checkpoint-step")->as<int>()
                  : n_iters,
        },
        .out_prefix = cmd.get_option("--out")->as<std::string>(),
    };

    if (cmd.get_option("--from-ckpt")->count() > 0)
    {
        config.from_checkpoint_path
            = cmd.get_option("--from-ckpt")->as<std::string>();
    }

    gelex::mcmc::ConfigValidator{config}.validate();
    return config;
}

}  // namespace cli
