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

#include <argparse.h>

#include "gelex/algo/mcmc/params.h"
#include "gelex/engine/mcmc.h"

namespace gelex::cli
{

auto make_mcmc_engine_config(argparse::ArgumentParser& cmd)
    -> mcmc::Engine::Config
{
    const auto n_iters = cmd.get<int>("--iters");
    mcmc::Engine::Config config{
        .bfile_prefix = cmd.get("--bfile"),
        .seed = cmd.get<int>("--seed"),
        .mcmc_params = gelex::mcmc::Params{
            .n_iters = n_iters,
            .n_burn_in = cmd.get<int>("--burn-in"),
            .n_thin = cmd.get<int>("--thin"),
            .checkpoint_step
            = cmd.is_used("--checkpoint-step")
                  ? cmd.get<int>("--checkpoint-step")
                  : n_iters,
        },
        .out_prefix = cmd.get("--out"),
    };

    if (cmd.is_used("--from-ckpt"))
    {
        config.from_checkpoint_path = cmd.get("--from-ckpt");
    }

    mcmc::ConfigValidator{config}.validate();
    return config;
}

}  // namespace gelex::cli
