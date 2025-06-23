#include "gelex/estimator/bayes/result.h"
#include <omp.h>
#include <cstddef>

#include "gelex/estimator/bayes/diagnostics.h"

namespace gelex
{
void MCMCResult::compute_summary_statistics(
    PosteriorSummary& summary,
    const arma::dcube& samples,
    double prob)
{
    if (samples.is_empty())
    {
        return;
    }
#pragma omp parallel for default(none) shared(samples, prob, summary)
    for (size_t r = 0; r < samples.n_rows; ++r)
    {
        arma::dvec flat_sample = arma::vectorise(samples.row(r));
        double mean = arma::mean(flat_sample);
        double median = arma::median(flat_sample);
        double stddev = arma::stddev(flat_sample);
        auto [hpdi_low, hpdi_high] = hpdi(flat_sample, prob);

        summary.mean.at(r) = mean;
        summary.median.at(r) = median;
        summary.stddev.at(r) = stddev;
        summary.hpdi_low.at(r) = hpdi_low;
        summary.hpdi_high.at(r) = hpdi_high;
    }
    summary.ess = effect_sample_size(samples, true);
    summary.rhat = split_gelman_rubin(samples);
}

MCMCResult::MCMCResult(const MCMCSamples& samples)
    : fixed(samples.fixed().n_rows), residual(1)
{
    for (size_t i = 0; i < samples.random().size(); ++i)
    {
        random.emplace_back(
            samples.random().coeffs[i].n_rows,
            samples.random().sigmas[i].n_rows);
    }

    for (size_t i = 0; i < samples.genetic().size(); ++i)
    {
        genetic.emplace_back(
            samples.genetic().coeffs[i].n_rows,
            samples.genetic().sigmas[i].n_rows,
            samples.genetic().pi[i].n_rows);
    }
}

void MCMCResult::compute_summary_statistics(
    const MCMCSamples& samples,
    double prob)
{
    compute_summary_statistics(fixed, samples.fixed(), prob);

    for (size_t i = 0; i < samples.random().size(); ++i)
    {
        compute_summary_statistics(
            random[i].coeff, samples.random().coeffs[i], prob);
        compute_summary_statistics(
            random[i].sigma, samples.random().sigmas[i], prob);
    }

    for (size_t i = 0; i < samples.genetic().size(); ++i)
    {
        compute_summary_statistics(
            genetic[i].pi, samples.genetic().pi[i], prob);
        compute_summary_statistics(
            genetic[i].genetic_var, samples.genetic().genetic_var[i], prob);
        compute_summary_statistics(
            genetic[i].heritability, samples.genetic().heritability[i], prob);
    }

    compute_summary_statistics(residual, samples.residual(), prob);
}

}  // namespace gelex
