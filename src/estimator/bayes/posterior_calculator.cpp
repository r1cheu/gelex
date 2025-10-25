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
 * For each parameter i:
 *   pve = (variance[i] * mean(beta_i)^2) / phenotype_var
 * Computes PVE using mean coefficients across all MCMC samples.
 *
 * @performance Computes PVE directly from mean coefficients, avoiding
 * iteration over individual MCMC samples for PVE calculation.
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
        double sum_beta = 0.0;

        // Compute mean beta across all MCMC samples
        for (const auto& chain_matrix : samples)
        {
            const auto row_view = chain_matrix.row(param_idx);
            sum_beta += row_view.sum();
        }

        // Compute mean beta and PVE
        const double mean_beta = sum_beta / total_draws;
        const double pve = var_ratio * mean_beta * mean_beta;

        summary.mean(param_idx) = pve;
        summary.stddev(param_idx) = 0.0;  // PVEse is no longer calculated
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

Eigen::MatrixXd compute_component_probs(
    const IntSamples& tracker_samples,
    Eigen::Index n_components)
{
    if (tracker_samples.empty() || tracker_samples[0].rows() == 0
        || n_components <= 0)
    {
        return Eigen::MatrixXd::Zero(0, 0);
    }

    const Eigen::Index n_snps = tracker_samples[0].rows();
    const auto n_chains = static_cast<Eigen::Index>(tracker_samples.size());
    const Eigen::Index n_draws = tracker_samples[0].cols();
    const double total_samples = static_cast<double>(n_chains * n_draws);

    Eigen::MatrixXd comp_probs = Eigen::MatrixXd::Zero(n_snps, n_components);

    // Count occurrences of each component for each SNP
    for (Eigen::Index chain = 0; chain < n_chains; ++chain)
    {
        for (int comp = 0; comp < n_components; ++comp)
        {
            // Count where tracker == comp for each SNP
            comp_probs.col(comp).array()
                += (tracker_samples[chain].array() == comp)
                       .cast<double>()
                       .rowwise()
                       .sum();
        }
    }

    // Normalize by total number of samples
    comp_probs /= total_samples;

    return comp_probs;
}

}  // namespace gelex::detail::PosteriorCalculator
