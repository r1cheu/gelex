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

#include "gelex/algo/infer/posterior_calculator.h"

#include <Eigen/Core>

namespace gelex::detail::PosteriorCalculator
{

PosteriorSummary compute_param_summary(
    const Eigen::Ref<const Eigen::MatrixXd>& samples,
    double prob)
{
    static_cast<void>(prob);

    if (samples.cols() == 0 || samples.rows() == 0)
    {
        return PosteriorSummary(0);
    }

    PosteriorSummary summary(get_n_params(samples));
    compute_mean_std(summary, samples);
    return summary;
}

PosteriorSummary compute_snp_summary(
    const Eigen::Ref<const Eigen::MatrixXd>& samples)
{
    if (samples.cols() == 0 || samples.rows() == 0)
    {
        return PosteriorSummary(0);
    }

    PosteriorSummary summary(get_n_params(samples));
    compute_mean_std(summary, samples);
    return summary;
}

/**
 * @brief Computes mean and standard deviation without intermediate memory
 * allocation.
 * @performance This is the most significant performance optimization. It avoids
 * flattening the samples by iterating directly over the original data,
 * dramatically reducing memory allocation and data copying.
 */
void compute_mean_std(
    PosteriorSummary& summary,
    const Eigen::Ref<const Eigen::MatrixXd>& samples)
{
    const Eigen::Index n_params = get_n_params(samples);
    if (n_params == 0)
    {
        return;
    }

    const auto total_draws = static_cast<double>(samples.cols());
    if (total_draws <= 1.0)
    {
        return;
    }

    const EigenThreadGuard guard;

#pragma omp parallel for default(none) \
    shared(summary, samples, n_params, total_draws)
    for (Eigen::Index param_idx = 0; param_idx < n_params; ++param_idx)
    {
        const auto row_view = samples.row(param_idx);
        const double sum = row_view.sum();
        const double sum_sq = row_view.array().square().sum();

        const double mean = sum / total_draws;
        const double variance
            = (sum_sq - total_draws * mean * mean) / (total_draws - 1.0);

        summary.mean(param_idx) = mean;
        summary.stddev(param_idx) = std::sqrt(std::max(0.0, variance));
    }
}

/**
 * @brief Computes PVE (proportion of variance explained) for each parameter.
 *
 * For each parameter i:
 *   pve = (variance[i] * mean(beta_i)^2) / phenotype_var
 * Computes PVE using mean coefficients across all MCMC samples.
 *
 * @performance Computes PVE directly from mean coefficients, avoiding
 * iteration over individual MCMC samples for PVE calculation.
 */
void compute_pve(
    PosteriorSummary& summary,
    const Eigen::Ref<const Eigen::MatrixXd>& samples,
    double phenotype_var)
{
    const Eigen::Index n_params = get_n_params(samples);
    if (n_params == 0 || phenotype_var <= 0.0)
    {
        return;
    }

    const auto total_draws = static_cast<double>(samples.cols());
    if (total_draws <= 1.0)
    {
        return;
    }

    const EigenThreadGuard guard;

#pragma omp parallel for default(none) \
    shared(summary, samples, n_params, total_draws, phenotype_var)
    for (Eigen::Index param_idx = 0; param_idx < n_params; ++param_idx)
    {
        const double mean_beta = samples.row(param_idx).sum() / total_draws;
        const double pve = mean_beta * mean_beta / phenotype_var;

        summary.mean(param_idx) = pve;
    }
}

Eigen::Index get_n_params(const Eigen::Ref<const Eigen::MatrixXd>& samples)
{
    return samples.rows();
}

Eigen::MatrixXd compute_component_probs(
    const Eigen::Ref<const Eigen::MatrixXi>& tracker_samples,
    Eigen::Index n_components)
{
    if (tracker_samples.size() == 0 || n_components <= 0)
    {
        return Eigen::MatrixXd::Zero(0, 0);
    }

    const Eigen::Index n_snps = tracker_samples.rows();
    const auto total_samples = static_cast<double>(tracker_samples.cols());

    Eigen::MatrixXd comp_probs = Eigen::MatrixXd::Zero(n_snps, n_components);

    for (int comp = 0; comp < n_components; ++comp)
    {
        comp_probs.col(comp).array()
            += (tracker_samples.array() == comp).cast<double>().rowwise().sum();
    }

    comp_probs /= total_samples;

    return comp_probs;
}

}  // namespace gelex::detail::PosteriorCalculator
