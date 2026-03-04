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

#include "gelex/algo/infer/reml.h"

#include <Eigen/Core>

#include "gelex/algo/infer/estimator.h"
#include "gelex/algo/numerics/optimizer_state.h"
#include "gelex/algo/numerics/variance_calculator.h"
#include "gelex/exception.h"
#include "gelex/infra/logger.h"
#include "gelex/infra/utils/formatter.h"
#include "gelex/model/freq/model.h"
#include "gelex/types/freq_effect.h"

namespace gelex
{

auto load_data_for_reml(
    const DataPipe::Config& config,
    DataPipeObserver observer) -> DataPipe
{
    DataPipe data_pipe(config, std::move(observer));
    data_pipe.load();
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
    logger->info(gelex::section("[Model Configuration]"));

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
