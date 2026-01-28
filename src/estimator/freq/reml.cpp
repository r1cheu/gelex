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

#include "gelex/estimator/freq/reml.h"
#include <Eigen/Core>

#include "../src/utils/formatter.h"
#include "gelex/estimator/freq/estimator.h"
#include "gelex/exception.h"
#include "gelex/logger.h"
#include "gelex/model/freq/model.h"
#include "gelex/optim/optimizer_state.h"
#include "gelex/optim/variance_calculator.h"
#include "types/freq_effect.h"

namespace gelex
{

auto load_data_for_reml(const DataPipe::Config& config) -> DataPipe
{
    auto logger = gelex::logging::get();
    DataPipe data_pipe(config);

    // Load phenotype first to get trait name
    auto p_stats = data_pipe.load_phenotypes();

    logger->info(gelex::section("Loading Data..."));
    logger->info(
        gelex::success(
            "Phenotypes : {:L} samples ({})",
            p_stats.samples_loaded,
            p_stats.trait_name));

    logger->info(
        gelex::success(
            "Genotypes  : {:L} samples", data_pipe.num_genotype_samples()));

    auto c_stats = data_pipe.load_covariates();
    if (c_stats.qcovar_loaded > 0 || c_stats.dcovar_loaded > 0)
    {
        logger->info(gelex::task("Covariates : "));
        if (c_stats.qcovar_loaded > 0)
        {
            logger->info(
                gelex::subtask(
                    "Quantitative : {} loaded ",
                    gelex::format_names(c_stats.q_names)));
        }
        if (c_stats.dcovar_loaded > 0)
        {
            logger->info(
                gelex::subtask(
                    "Discrete     : {} loaded ",
                    gelex::format_names(c_stats.d_names)));
        }
    }

    // Load GRM(s)
    logger->info(gelex::success("GRM : "));
    auto grm_stats = data_pipe.load_grms();
    std::string grm_str;

    for (const auto& grm_stat : grm_stats)
    {
        switch (grm_stat.type)
        {
            case freq::GrmType::A:
                logger->info(
                    gelex::subtask(
                        "Additive     : {:L} samples",
                        grm_stat.samples_in_file));
                grm_str += "Additive, ";
                break;
            case freq::GrmType::D:
                logger->info(
                    gelex::subtask(
                        "Dominance    : {:L} samples",
                        grm_stat.samples_in_file));
                grm_str += "Dominance, ";
                break;
            default:
                logger->info(
                    gelex::subtask(
                        "Unknown      : {:L} samples",
                        grm_stat.samples_in_file));
                grm_str += "Unknown, ";
                break;
        }
    }

    logger->info("");
    logger->info(gelex::section("Pre-processing..."));
    auto i_stats = data_pipe.intersect_samples();
    logger->info(gelex::task("Sample Intersection:"));
    logger->info(gelex::subtask("Common   : {:L}", i_stats.common_samples));
    logger->info(gelex::subtask("Excluded : {:L}", i_stats.excluded_samples));

    if (i_stats.common_samples == 0)
    {
        throw gelex::InvalidInputException(
            "No common samples found between phenotype, covariates, GRM");
    }

    data_pipe.finalize();
    return data_pipe;
}

auto reml(
    const DataPipe::Config& config,
    size_t max_iter,
    double tol,
    bool em_init,
    bool verbose) -> std::
    tuple<std::shared_ptr<SampleManager>, Eigen::MatrixXd, Eigen::VectorXd>
{
    auto logger = gelex::logging::get();
    auto data_pipe = load_data_for_reml(config);

    logger->info("");
    logger->info(gelex::section("Model Configuration..."));

    gelex::FreqModel model(data_pipe);
    gelex::FreqState state(model);

    logger->info(gelex::task("Design:"));
    logger->info(
        gelex::subtask("Fixed Effects   : {}", model.fixed().X.cols()));

    std::string grm_str;
    for (const auto& g : model.genetic())
    {
        if (g.type == freq::GrmType::A)
        {
            grm_str += "Additive, ";
        }
        else if (g.type == freq::GrmType::D)
        {
            grm_str += "Dominance, ";
        }
        else
        {
            grm_str += "Unknown, ";
        }
    }
    if (!grm_str.empty())
    {
        grm_str = grm_str.substr(0, grm_str.length() - 2);
    }

    logger->info(
        gelex::subtask(
            "Genetic Effects : {} ({})", model.genetic().size(), grm_str));

    logger->info(gelex::task("Optimizer (AI):"));
    logger->info(gelex::subtask("Tolerance : {:.1e}", tol));
    logger->info(gelex::subtask("Max Iter  : {}", max_iter));

    logger->info("");
    logger->info(gelex::section("Fitting Null Model..."));

    gelex::Estimator estimator(max_iter, tol);
    Eigen::MatrixXd v = estimator.fit(model, state, em_init, verbose);
    Eigen::VectorXd v_inv_residual
        = v * (model.phenotype() - model.fixed().X * state.fixed().coeff);

    if (!estimator.is_converged())
    {
        logger->warn("REML did not converge, results may be unreliable");
    }
    return {
        data_pipe.sample_manager(), std::move(v), std::move(v_inv_residual)};
}

}  // namespace gelex
