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

#include "vi_config.h"

#include <argparse.h>

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
        .dominance = bayes::DominancePolicy::symmetric,
        .estimate_pi = false,
    };

    if (base != BayesBase::RR)
    {
        throw gelex::GelexException(
            fmt::format("CAVI only supports BayesRR, got: {}", method));
    }
    if (!bayes::is_valid_method(method))
    {
        throw gelex::GelexException(
            fmt::format("invalid method combination: {}", method));
    }
    return method;
}

}  // namespace

auto make_vi_config(argparse::ArgumentParser& cmd) -> vi::Engine::Config
{
    auto method = parse_method(cmd);
    return vi::Engine::Config{
        .bfile_prefix = cmd.get("--bfile"),
        .method = method,
        .params = vi::Params{
            .max_iters = cmd.get<int>("--max-iters"),
            .tol = cmd.get<double>("--tol"),
        },
        .out_prefix = cmd.get("--out"),
    };
}

}  // namespace gelex::cli
