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

#include "cli/cli_helper.h"
#include "gelex/data/genotype/bed_path.h"
#include "gelex/exception.h"

namespace
{

auto has_dominance(gelex::BayesAlphabet type) -> bool
{
    switch (type)
    {
        case gelex::BayesAlphabet::Bd:
        case gelex::BayesAlphabet::Bdpi:
        case gelex::BayesAlphabet::Cd:
        case gelex::BayesAlphabet::Cdpi:
        case gelex::BayesAlphabet::Rd:
        case gelex::BayesAlphabet::Ad:
        case gelex::BayesAlphabet::RRd:
            return true;
        default:
            return false;
    }
}

}  // namespace

auto FitConfig::make(argparse::ArgumentParser& cmd) -> FitConfig
{
    auto method_name = cmd.get("-m");
    auto method = gelex::get_bayesalphabet(method_name)
                      .value_or(gelex::BayesAlphabet::RR);
    auto use_dominance = has_dominance(method);

    FitConfig config{
        .method_name = method_name,
        .method = method,
        .use_dominance = use_dominance,
        .threads = cmd.get<int>("--threads"),
        .seed = cmd.get<int>("--seed"),
        .out_prefix = cmd.get("--out"),
        .mcmc_params = gelex::MCMCParams(
            cmd.get<int>("--iters"),
            cmd.get<int>("--burnin"),
            cmd.get<int>("--thin")),
        .phenotype_path = cmd.get("pheno"),
        .phenotype_column = cmd.get<int>("--pheno-col"),
        .bed_path = gelex::format_bed_path(cmd.get("bfile")),
        .use_mmap = cmd.get<bool>("--mmap"),
        .chunk_size = cmd.get<int>("--chunk-size"),
        .qcovar_path = cmd.get("--qcovar"),
        .dcovar_path = cmd.get("--dcovar"),
        .genotype_method = gelex::cli::parse_genotype_process_method(
            cmd.get<std::string>("--geno-method"))};

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

    return config;
}

auto FitConfig::validate() const -> void
{
    if (chunk_size <= 0)
    {
        throw gelex::InvalidInputException("chunk_size must be > 0");
    }
    if (threads != -1 && threads <= 0)
    {
        throw gelex::InvalidInputException("threads must be -1 or > 0");
    }
    if (mcmc_params.n_burnin >= mcmc_params.n_iters)
    {
        throw gelex::InvalidInputException(
            "n_burnin must be less than n_iters");
    }
}
