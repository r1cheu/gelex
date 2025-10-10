#include "posterior_calculator.h"

#include <Eigen/Core>

#include "estimator/bayes/diagnostics.h"
#include "gelex/estimator/bayes/samples.h"

namespace gelex::detail::PosteriorCalculator
{

PosteriorSummary compute_param_summary(const Samples& samples, double prob)
{
    if (samples.empty() || samples[0].rows() == 0)
    {
        return PosteriorSummary(0);
    }

    PosteriorSummary summary(get_n_params(samples));
    compute_mean_std(summary, samples);
    compute_hpdi(summary, samples, prob);
    compute_ess(summary, samples);
    compute_rhat(summary, samples);
    return summary;
}

PosteriorSummary compute_snp_summary(const Samples& samples)
{
    if (samples.empty() || samples[0].rows() == 0)
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
void compute_mean_std(PosteriorSummary& summary, const Samples& samples)
{
    const Eigen::Index n_params = get_n_params(samples);
    if (n_params == 0)
    {
        return;
    }

    const auto total_draws
        = static_cast<double>(get_n_chains(samples) * get_n_draws(samples));
    if (total_draws <= 1.0)
    {
        return;
    }

    const EigenThreadGuard guard;

#pragma omp parallel for default(none) \
    shared(summary, samples, n_params, total_draws)
    for (Eigen::Index param_idx = 0; param_idx < n_params; ++param_idx)
    {
        double sum = 0.0;
        double sum_sq = 0.0;

        for (const auto& chain_matrix : samples)
        {
            const auto row_view = chain_matrix.row(param_idx);
            sum += row_view.sum();
            sum_sq += row_view.array().square().sum();
        }

        const double mean = sum / total_draws;
        const double variance
            = (sum_sq - total_draws * mean * mean) / (total_draws - 1.0);

        summary.mean(param_idx) = mean;
        summary.stddev(param_idx) = std::sqrt(std::max(0.0, variance));
    }
}

void compute_hpdi(
    PosteriorSummary& summary,
    const Samples& samples,
    double prob)
{
    if (samples.empty() || samples[0].rows() == 0)
    {
        return;
    }

    const Eigen::Index n_params = get_n_params(samples);
    const Eigen::Index n_chains = get_n_chains(samples);
    const Eigen::Index n_draws = get_n_draws(samples);

    const EigenThreadGuard guard;

#pragma omp parallel for default(none) \
    shared(samples, prob, summary, n_params, n_chains, n_draws)
    for (Eigen::Index param_idx = 0; param_idx < n_params; ++param_idx)
    {
        Eigen::VectorXd flat_sample = flatten_samples(samples, param_idx);
        auto [hpdi_low_val, hpdi_high_val] = hpdi(flat_sample, prob);

        summary.hpdi_low(param_idx) = hpdi_low_val;
        summary.hpdi_high(param_idx) = hpdi_high_val;
    }
}

void compute_ess(PosteriorSummary& summary, const Samples& samples)
{
    if (samples.empty() || samples[0].rows() == 0)
    {
        return;
    }

    summary.ess = effect_sample_size(samples, true);
}

void compute_rhat(PosteriorSummary& summary, const Samples& samples)
{
    if (samples.empty() || samples[0].rows() == 0)
    {
        return;
    }

    summary.rhat = split_gelman_rubin(samples);
}

/**
 * @brief Computes PVE (proportion of variance explained) for each parameter.
 *
 * For each parameter i and each MCMC sample:
 *   pve_sample = (variance[i] * beta_sample[i]^2) / phenotype_var
 * Then computes mean and std of PVE across all samples.
 *
 * @performance Precalculates variance ratios and uses direct iteration over
 * chains and draws without flattening, parallelizes across parameters with
 * OpenMP.
 */
void compute_pve(
    PosteriorSummary& summary,
    const Samples& samples,
    const Eigen::VectorXd& variances,
    double phenotype_var)
{
    const Eigen::Index n_params = get_n_params(samples);
    if (n_params == 0 || phenotype_var <= 0.0)
    {
        return;
    }

    if (variances.size() != n_params)
    {
        return;
    }

    const auto total_draws
        = static_cast<double>(get_n_chains(samples) * get_n_draws(samples));
    if (total_draws <= 1.0)
    {
        return;
    }

    // Precalculate variance ratios
    Eigen::VectorXd var_ratios = variances / phenotype_var;

    const EigenThreadGuard guard;

#pragma omp parallel for default(none) \
    shared(summary, samples, var_ratios, n_params, total_draws)
    for (Eigen::Index param_idx = 0; param_idx < n_params; ++param_idx)
    {
        const double var_ratio = var_ratios(param_idx);
        double sum_pve = 0.0;
        double sum_pve_sq = 0.0;

        // Compute PVE for each MCMC sample
        for (const auto& chain_matrix : samples)
        {
            const auto row_view = chain_matrix.row(param_idx);
            for (Eigen::Index draw = 0; draw < row_view.size(); ++draw)
            {
                const double beta = row_view(draw);
                const double pve_sample = var_ratio * beta * beta;
                sum_pve += pve_sample;
                sum_pve_sq += pve_sample * pve_sample;
            }
        }

        // Compute mean and std of PVE
        const double mean_pve = sum_pve / total_draws;
        const double variance_pve
            = (sum_pve_sq - total_draws * mean_pve * mean_pve)
              / (total_draws - 1.0);

        summary.mean(param_idx) = mean_pve;
        summary.stddev(param_idx) = std::sqrt(std::max(0.0, variance_pve));
    }
}

Eigen::VectorXd flatten_samples(
    const Samples& samples,
    Eigen::Index param_index)
{
    const Eigen::Index n_chains = get_n_chains(samples);
    const Eigen::Index n_draws = get_n_draws(samples);

    Eigen::VectorXd flat_sample(n_chains * n_draws);
    for (Eigen::Index chain = 0; chain < n_chains; ++chain)
    {
        flat_sample.segment(chain * n_draws, n_draws)
            = samples[chain].row(param_index).transpose();
    }
    return flat_sample;
}

Eigen::Index get_n_params(const Samples& samples)
{
    return samples.empty() ? 0 : samples[0].rows();
}

Eigen::Index get_n_chains(const Samples& samples)
{
    return static_cast<Eigen::Index>(samples.size());
}

Eigen::Index get_n_draws(const Samples& samples)
{
    return samples.empty() ? 0 : samples[0].cols();
}

}  // namespace gelex::detail::PosteriorCalculator
