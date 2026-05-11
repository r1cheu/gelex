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

#include <vector>

#include <argparse.h>

#include "gelex/exception.h"
#include "gelex/types/genetic_effect_type.h"

namespace gelex::cli
{

namespace
{

auto parse_mode_flag(std::string_view sv) -> std::vector<GeneticMode>
{
    if (sv == "A")
    {
        return {GeneticMode::A};
    }
    if (sv == "D")
    {
        return {GeneticMode::D};
    }
    if (sv == "AD")
    {
        return {GeneticMode::A, GeneticMode::D};
    }
    throw GelexException(fmt::format("invalid --mode: {}", sv));
}

auto parse_method(argparse::ArgumentParser& cmd)
    -> std::pair<bayes::BayesConfig, std::vector<GeneticMode>>
{
    auto base
        = gelex::get_bayes_base(cmd.get("-m")).value_or(gelex::BayesBase::RR);
    auto requested = parse_mode_flag(cmd.get("--mode"));

    bayes::BayesConfig cfg{
        .base = base,
        .dominance = bayes::DominancePolicy::symmetric,
        .estimate_pi = false,
    };

    if (base != BayesBase::RR)
    {
        throw gelex::GelexException(
            fmt::format("CAVI only supports BayesRR, got: {}", cfg));
    }
    if (!bayes::is_valid_method(cfg, requested))
    {
        throw gelex::GelexException(
            fmt::format("invalid method combination: {}", cfg));
    }
    return {cfg, std::move(requested)};
}

}  // namespace

auto make_vi_config(argparse::ArgumentParser& cmd) -> vi::Engine::Config
{
    auto [method, requested] = parse_method(cmd);
    return vi::Engine::Config{
        .bfile_prefix = cmd.get("--bfile"),
        .method = method,
        .requested_effects = std::move(requested),
        .params = vi::Params{
            .max_iters = cmd.get<int>("--max-iters"),
            .tol = cmd.get<double>("--tol"),
        },
        .out_prefix = cmd.get("--out"),
    };
}

}  // namespace gelex::cli
