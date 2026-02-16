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

#ifndef GELEX_ESTIMATOR_BAYES_POSTERIOR_CALCULATOR_H_
#define GELEX_ESTIMATOR_BAYES_POSTERIOR_CALCULATOR_H_

#include <Eigen/Core>

#include "gelex/types/mcmc_results.h"
#include "gelex/types/mcmc_samples.h"

namespace gelex::detail
{

class EigenThreadGuard
{
   public:
    EigenThreadGuard() : old_thread_count_(Eigen::nbThreads())
    {
        Eigen::setNbThreads(1);
    }

    ~EigenThreadGuard() { Eigen::setNbThreads(old_thread_count_); }

    EigenThreadGuard(const EigenThreadGuard&) = delete;
    EigenThreadGuard& operator=(const EigenThreadGuard&) = delete;
    EigenThreadGuard(EigenThreadGuard&&) = delete;
    EigenThreadGuard& operator=(EigenThreadGuard&&) = delete;

   private:
    int old_thread_count_;
};

namespace PosteriorCalculator
{

PosteriorSummary compute_param_summary(const Samples& samples, double prob);

PosteriorSummary compute_snp_summary(const Samples& samples);

void compute_mean_std(PosteriorSummary& summary, const Samples& samples);

void compute_pve(
    PosteriorSummary& summary,
    const Samples& samples,
    double phenotype_var);

Eigen::VectorXd flatten_samples(
    const Samples& samples,
    Eigen::Index param_index);

Eigen::Index get_n_params(const Samples& samples);

Eigen::Index get_n_chains(const Samples& samples);

Eigen::Index get_n_draws(const Samples& samples);

Eigen::MatrixXd compute_component_probs(
    const IntSamples& tracker_samples,
    Eigen::Index n_components);

}  // namespace PosteriorCalculator

}  // namespace gelex::detail

#endif  // GELEX_ESTIMATOR_BAYES_POSTERIOR_CALCULATOR_H_
